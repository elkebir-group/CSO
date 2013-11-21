/* 
 * main.cpp
 *
 *  Created on: 26-apr-2011
 *      Author: M. El-Kebir
 */

#include <iostream>
#include <fstream>
#include <string>
#include <lemon/arg_parser.h>
#include <lemon/time_measure.h>
#include "common/data.h"
#include "ilpsolver.h"
#include "ilpsolverext.h"
#include "analysis/lowerbounds.h"
#include <sstream>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>


void handler(int sig)
{
  void *array[30];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 30);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, 2);
  exit(1);
}

IlpSolver::SolverStatus solve(const Data* pData,
                              const IlpSolver::Options& options,
                              bool feasiblity,
                              int timeLimit,
                              CrossingSchedule*& pBestSchedule,
                              bool pareto,
                              bool correctProb)
{
  IlpSolver* pSolver1 = correctProb ? new IlpSolverExt(pData, options) :
                                      new IlpSolver(pData, options);
  IlpSolver* pSolver2 = correctProb ? new IlpSolverExt(pData, options) :
                                      new IlpSolver(pData, options);
  bool solve1 = false;
  bool solve2 = false;
  IlpSolver::SolverStatus stat1, stat2, finalStat;

  if (timeLimit < 0)
  {
    finalStat = IlpSolver::CSO_SOLVER_INFEASIBLE;

    pSolver1->init(true);
    stat1 = pSolver1->solve(feasiblity, -1);
    solve1 =  stat1 == IlpSolver::CSO_SOLVER_OPTIMAL;

    if (!pData->isIdeotypeHomozygous())
    {
      pSolver2->init(false);
      stat2 = pSolver2->solve(feasiblity, -1);
      solve2 = stat2 == IlpSolver::CSO_SOLVER_OPTIMAL;
    }
  }
  else
  {
    finalStat = IlpSolver::CSO_SOLVER_INFEASIBLE;

    pSolver1->init(true);
    stat1 = pSolver1->solve(feasiblity, timeLimit);
    solve1 = stat1 != IlpSolver::CSO_SOLVER_INFEASIBLE
        && stat1 != IlpSolver::CSO_SOLVER_TIME_LIMIT_INFEASIBLE;

    if (!pData->isIdeotypeHomozygous())
    {
      pSolver2->init(false);
      stat2 = pSolver2->solve(feasiblity, timeLimit);
      solve2 = stat2 != IlpSolver::CSO_SOLVER_INFEASIBLE
          && stat2 != IlpSolver::CSO_SOLVER_TIME_LIMIT_INFEASIBLE;
    }
  }

  if (solve1 || solve2)
  {
    if (solve1 && solve2)
    {
      if (pSolver1->getObjectiveValue() < pSolver2->getObjectiveValue())
      {
        finalStat = stat1;
        if (pBestSchedule && pBestSchedule->getCost() > pSolver1->getObjectiveValue())
        {
          delete pBestSchedule;
          pBestSchedule = new CrossingSchedule(*pSolver1);
        }
        else if (!pBestSchedule || pareto)
        {
          pBestSchedule = new CrossingSchedule(*pSolver1);
        }
      }
      else
      {
        finalStat = stat2;
        if (pBestSchedule && pBestSchedule->getCost() > pSolver2->getObjectiveValue())
        {
          delete pBestSchedule;
          pBestSchedule = new CrossingSchedule(*pSolver2);
        }
        else if (!pBestSchedule || pareto)
        {
          pBestSchedule = new CrossingSchedule(*pSolver2);
        }
      }
    }
    else if (solve1)
    {
      finalStat = stat1;
      if (pBestSchedule && pBestSchedule->getCost() > pSolver1->getObjectiveValue())
      {
        delete pBestSchedule;
        pBestSchedule = new CrossingSchedule(*pSolver1);
      }
      else if (!pBestSchedule || pareto)
      {
        pBestSchedule = new CrossingSchedule(*pSolver1);
      }
    }
    else
    {
      finalStat = stat2;
      if (pBestSchedule && pBestSchedule->getCost() > pSolver2->getObjectiveValue())
      {
        delete pBestSchedule;
        pBestSchedule = new CrossingSchedule(*pSolver2);
      }
      else if (!pBestSchedule || pareto)
      {
        pBestSchedule = new CrossingSchedule(*pSolver2);
      }
    }
  }

  delete pSolver1;
  delete pSolver2;

  return finalStat;
}

int main(int argc, char** argv)
{
#ifdef NDEBUG
  // crash gracefully in release mode
  signal(SIGSEGV, handler);   // install our handler
  signal(SIGABRT, handler);   // install our handler
#endif

  IlpSolver::Options options;
  bool printRM = false;
  std::string outputFileName;
  bool automatic = false;
  int timeLimit = -1;
  bool saveIntermediate = false;
  bool pareto = false;
  bool printProb = false;
  bool correctProb = false;
  int nInternalNodes = 0;
  int nGen = 0;

  lemon::ArgParser ap(argc, argv);
  ap.refOption("v", "Verbose CPLEX output", 
        options._verbose, false)
    .refOption("pareto", "Pareto", pareto, false)
    //.refOption("t", "Enforce crossing schedules to be trees",
    //    options._tree, false)
    .refOption("cc", "Correct probability function", correctProb, false)
    .refOption("ns", "Cuts: no selfing",
        options._noSelfing, false)
    .refOption("c", "Cuts: useful crossovers only",
        options._usefulCross, false)
    .refOption("n", "Cuts: based on bound Nmax",
        options._boundNmax, false)
    .refOption("b", "Number of internal nodes",
        nInternalNodes, false)
    .refOption("br", "Generation",
        nGen, false)
    .refOption("ub", "Upper bound on the objective value",
        options._upperBoundObj, false)
    .refOption("ub", "Upper bound on the total population size",
        options._upperBoundPop, false)
    .refOption("servin", "Servin instances",
        options._servin, false)
    .refOption("intermediate", "Save intermediate results (use in conjuction with -a)",
        saveIntermediate, false)
    .refOption("timeLimit", "Time limit in seconds (use in conjuction with -a)",
        timeLimit, false)
    .refOption("g", "Number of segments used to approximate population size",
        options._nof_segs, false)
    .refOption("print-prob", "Print prob on edges",
        printProb, false)
    .refOption("print-RM", "Print RM", printRM, false)
    .refOption("o", "Output file name",
        outputFileName, false)
    .refOption("a", "Automatic mode", automatic, false);
  ap.parse();

  if (nInternalNodes < 0)
  {
    std::cerr << "Number of internal nodes must be a nonzero integer" << std::endl;
    return 1;
  }
  options._bound = static_cast<size_t>(nInternalNodes);

  if (nGen < 0)
  {
    std::cerr << "Number of generations must be a nonzero integer" << std::endl;
    return 1;
  }
  options._fixedGen = static_cast<size_t>(nGen);

  if (ap.files().size() == 0)
  {
    std::cerr << "Missing input file" << std::endl;
    return 1;
  }

  if (!ap.given("a") && !ap.given("b") && !ap.given("br"))
  {
    std::cerr << "Specify either -a or (-b and -br)" << std::endl;
    return 1;
  }

  const std::string& inputFileName = ap.files()[0];
  Data* pData = Data::create(inputFileName.c_str());
  if (!pData) return 1;

  if (printRM)
  {
    pData->printRM();
    pData->printLogRM();
    return 0;
  }
 
  if (options._servin)
  {
    options._bound = pData->getNumberOfLoci()+1;
  }

  bool overallOptimal = false;
  bool timeLimitHit = false;
  CrossingSchedule* pBestSchedule = NULL;
  lemon::Timer t;
  if (!automatic)
  {
    pData->updateGamma(options._bound);
    solve(pData, options, false, timeLimit, pBestSchedule, pareto, correctProb);
  }
  else
  {
    LowerBound LB(pData, false);
    bool foundFeasible = false;
    options._upperBoundPop = std::numeric_limits<double>::max();
    for (int b = std::max(LB.getCrossLB(), options._bound);; b++)
    {
      const int lbr = 1+ceil(log(b)/log(2));
      if (lbr > pData->getGenMax())
        break;
      pData->updateGamma(b);

      if (!pareto && pData->getCost(LB.getPopLB(), lbr, b) > options._upperBoundObj)
      {
        overallOptimal = true;
        std::cerr << "No need to look any further (" << b << "," << lbr << ")"
          << " as " << pData->getCost(LB.getPopLB(), lbr, b) 
          << " > " << options._upperBoundObj << std::endl;

        break;
      }

      // This doesn't not work, as _upperBoundPop may be computed using a br > lbr
      //if (LB.getPopLB() > options._upperBoundPop)
      //{
      //  overallOptimal = true;
      //  std::cerr << "No need to look any further (" << b << "," << lbr << ")"
      //    << " as " << LB.getPopLB() 
      //    << " > " << options._upperBoundPop << std::endl;

      //  break;
      //}

      for (int br = lbr; br <= std::min(b, pData->getGenMax()); br++)
      {
        if (!pareto)
        {
          options._upperBoundPop = (options._upperBoundObj
            - b * pData->getCostNode()
            - br * pData->getCostGen()) / pData->getCostCrossover();
        }

        std::cerr << "Solving (" << b << "," << br  << ")... " 
          << "pop <= " << options._upperBoundPop << " " << std::flush;

        options._bound = b;
        options._fixedGen = br;
        lemon::Timer t2;
        IlpSolver::SolverStatus stat = solve(pData,
                                             options,
                                             !foundFeasible,
                                             foundFeasible ? timeLimit : -1,
                                             pBestSchedule,
                                             pareto,
                                             correctProb);

        if (stat == IlpSolver::CSO_SOLVER_TIME_LIMIT_FEASIBLE)
        {
          stat = solve(pData,
                       options,
                       false,
                       -1,
                       pBestSchedule,
                       pareto,
                       correctProb);
          assert(stat == IlpSolver::CSO_SOLVER_OPTIMAL);
        }

        if (stat == IlpSolver::CSO_SOLVER_OPTIMAL)
        {
          if (!pareto)
            options._upperBoundObj = pBestSchedule->getCost();

          foundFeasible = true;

          std::cerr << "Done... Obj value: " << options._upperBoundObj << std::endl;

          fprintf(stderr, "\"%s\",%.3f,%lu,%lu,%.2f,%.2f,\"%s\"\n", 
            inputFileName.c_str(), t2.realTime(),
            pBestSchedule->getGen(), pBestSchedule->getCross(), pBestSchedule->getPop(),
            pBestSchedule->getCost(), "OPT");

          if (saveIntermediate && ap.given("o"))
          {
            std::stringstream ss;
            ss << outputFileName.substr(0, outputFileName.length() - 4) 
              << "_" << b << "_" << br << ".lgf";
            std::ofstream os(ss.str().c_str());
            if (os)
              pBestSchedule->saveDAG(os);
            os.flush();
            os.close();
          }
        }
        else
        {
          timeLimitHit |= stat == IlpSolver::CSO_SOLVER_TIME_LIMIT_INFEASIBLE;
          std::cerr << "Infeasible" << std::endl;
          fprintf(stderr, "\"%s\",%.3f,-,-,-,-,\"%s\"\n", inputFileName.c_str(), t2.realTime(),
              stat == IlpSolver::CSO_SOLVER_INFEASIBLE ? "INFEASIBLE" : "TIME LIMIT");
        }
      }
    }
  }
  t.stop();

  if (pBestSchedule)
  {
    if (ap.given("o"))
    {
      std::ofstream os(outputFileName.c_str());
      if (os)
        pBestSchedule->saveDAG(os);
    }
    else
    {
      pBestSchedule->printDAG(std::cout, printProb);
    }

    fprintf(stderr, "\"%s\",%.3f,%lu,%lu,%.2f,%.2f,\"%s\",", 
      inputFileName.c_str(), t.realTime(),
      pBestSchedule->getGen(), pBestSchedule->getCross(), pBestSchedule->getPop(),
      pBestSchedule->getCost(), !timeLimitHit && overallOptimal ? "OPT" : "?");
    pData->updateGamma(pBestSchedule->getCross());
    pBestSchedule->recomputePop();
    fprintf(stderr, "%.3f,%.3f\n", pBestSchedule->getPop(), pBestSchedule->getCost());

    delete pData;
    delete pBestSchedule;
    return 0;
  }
  else
  {
    fprintf(stderr, "\"%s\",%.3f,-,-,-,-,-\n", inputFileName.c_str(), t.realTime());

    delete pData;
    return 1;
  }
}
