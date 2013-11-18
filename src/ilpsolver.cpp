/*
 * ilpsolver.cpp
 *
 *  Created on: 26-apr-2011
 *      Author: M. El-Kebir
 */

#include "ilpsolver.h"
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include "analysis/lowerbounds.h"
#include "ilpsolverheuristic.h"

IlpSolver::IlpSolver(const Data* pData, const Options& options)
  : CrossingSchedule(pData)
  , _options(options)
  , _tol(1e-8)
  , _idx(_G)
  , _env()
  , _model(_env)
  , _allVar(_env)
  , _obj(_env)
  , _pCplex(NULL)
  , _breakpoint()
  , _N()
  , _cs1(_env, pData->getNumberOfLoci())
  , _cs2(_env, pData->getNumberOfLoci())
  , _c(_env, 2 * pData->getParents().size())
  , _homozygous(_env, pData->getParents().size())
  , _x(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _g(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _y(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _a(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _at(_env, pData->isIdeotypeHomozygous() ? options._bound-1 : options._bound)
  , _h(_env, pData->isIdeotypeHomozygous() ? options._bound-1 : options._bound)
  , _hx(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _d(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _p(_env, options._bound, -IloInfinity, 0.0, IloNumVar::Float)
  , _lambda(_env, options._bound)
  , _r(_env, options._bound, 1, options._fixedGen)
  , _b(_env, pData->isIdeotypeHomozygous() ? options._bound-1 : options._bound)
  , _bx(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _bxz(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _z(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _uc(_env, options._bound)
  , _ucx(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
{
}

IlpSolver::~IlpSolver()
{
  delete _pCplex;
  _obj.end();
  _cs1.end();
  _cs2.end();
  _c.end();
  _homozygous.end();
  _x.end();
  _g.end();
  _y.end();
  _a.end();
  _at.end();
  _h.end();
  _hx.end();
  _d.end();
  _p.end();
  _lambda.end();
  _r.end();
  _b.end();
  _bx.end();
  _bxz.end();
  _z.end();
  _uc.end();
  _ucx.end();
}

IlpSolver::SolverStatus IlpSolver::solve(bool feasibility, int timeLimit)
{
  delete _pCplex;
  _pCplex = new IloCplex(_model);

  if (!_options._verbose)
  {
    _pCplex->setOut(_env.getNullStream());
    _pCplex->setWarning(_env.getNullStream());
    _pCplex->setError(_env.getNullStream());
  }
  else
  {
    _pCplex->setOut(std::cout);
    _pCplex->setWarning(std::cout);
    _pCplex->setError(std::cout);
  }

  //_pCplex->setParam(IloCplex::AggFill, 0);
  //_pCplex->setParam(IloCplex::PreInd, 0);
  //_pCplex->setParam(IloCplex::RelaxPreInd, 0);
  //_pCplex->setParam(IloCplex::PreslvNd, -1);
  //_pCplex->setParam(IloCplex::RepeatPresolve, 0);

  //_pCplex->use(new IlpSolverHeuristic(_env,
  //                                    static_cast<int>(_pData->getParents().size()),
  //                                    _pData->getNumberOfLoci(),
  //                                    _options._fixedGen,
  //                                    _options._bound,
  //                                    _c,
  //                                    _obj,
  //                                    _x,
  //                                    _p,
  //                                    _allVar));

  if (timeLimit>0)
    _pCplex->setParam(IloCplex::TiLim, timeLimit);
  
  if (feasibility)
    _pCplex->setParam(IloCplex::MIPEmphasis, CPX_MIPEMPHASIS_FEASIBILITY);

  _pCplex->setParam(IloCplex::Threads, 1);
  //_pCplex->setParam(IloCplex::AggFill, 0);
  //_pCplex->setParam(IloCplex::PreInd, 0);
  //_pCplex->setParam(IloCplex::RelaxPreInd, 0);
  //_pCplex->setParam(IloCplex::PreslvNd, -1);
  //_pCplex->setParam(IloCplex::RepeatPresolve, 0);

  //_model.add(_x[0][3] == 1);
  //_model.add(_x[1][1] == 1);
  //_model.add(_x[2][4] == 1);
  //_model.add(_x[3][2] == 1);
  //_model.add(_a[0][3] == 0);
  //_model.add(_a[0][2] == 0);
  //_model.add(_a[0][1] == 1);
  //_model.add(_a[0][0] == 0);
  //_model.add(_a[1][3] == 1);
  //_model.add(_a[1][2] == 0);
  //_model.add(_a[1][1] == 0);
  //_model.add(_a[1][0] == 1);
  //_model.add(_a[2][3] == 1);
  //_model.add(_a[2][2] == 0);
  //_model.add(_a[2][1] == 1);
  //_model.add(_a[2][0] == 1);
  //_model.add(_d[2][2] == 1);
  //char buf[1024];
  //sprintf(buf, "letssee-c%lu-g%lu.lp", _options._bound, _options._fixedGen);
  //_pCplex->exportModel(buf);
  //std::cout << "// crs: " << _options._bound << std::endl;
  //std::cout << "// gen: " << _options._fixedGen << std::endl;
  //std::cout << "// Number of cols: " << _pCplex->getNcols() << " " << _allVar.getSize() << std::endl;
  //std::cout << "// Number of rows: " << _pCplex->getNrows() << std::endl;

  bool res = false;

  try
  {
    res = _pCplex->solve();
  }
  catch (IloCplex::Exception& e)
  {
    std::cerr << e.getStatus() << std::endl;
    exit(1);
  }

  if (!res)
  {
    //_env.error() << "Failed to optimize LP" << std::endl;
    return _pCplex->getCplexStatus() == CPX_STAT_ABORT_TIME_LIM ? 
      CSO_SOLVER_TIME_LIMIT_INFEASIBLE : CSO_SOLVER_INFEASIBLE;
  }
  else
  {
    if (_options._verbose)
    {
      std::cout << "// Objective value: " << _pCplex->getObjValue() << std::endl;
      //printInnerNodes();
      //std::cout << "//" << std::endl;
      //printA();
      //std::cout << "//" << std::endl;
      //printX();
      //std::cout << "//" << std::endl;
      //printY();
      //std::cout << "//" << std::endl;
      //printAt();
      //std::cout << "//" << std::endl;
      //printP();
      //std::cout << "//" << std::endl;
      //printLambda();
      //std::cout << "//" << std::endl;
      //printB();
      //std::cout << "//" << std::endl;
      //printZ();
      //std::cout << "//" << std::endl;
      //printG();
      //std::cout << "//" << std::endl;
      //printGG();
      //std::cout << "//" << std::endl;
      //printE();
      //std::cout << "//" << std::endl;
      //printHX();
      //std::cout << "//" << std::endl;
      //printHXE();
    }
    constructDAG();

    if (_pCplex->getCplexStatus() == CPX_STAT_ABORT_TIME_LIM)
    {
      return CSO_SOLVER_TIME_LIMIT_FEASIBLE;
    }

    return CSO_SOLVER_OPTIMAL;
  }
}

void IlpSolver::init(bool swapIdeotype)
{
  initSegments();
  initParentalGenotypes();
  initChromosomesFromGenotypes();
  initAllelesFromChromosomes();
  initIdeotype(swapIdeotype);

  initGen();
  initPop();
  initObj();

  _model.add(IloMinimize(_env, _obj));

  //if (_options._usefulCross)
  //  initUsefulCross();

  if (_options._noSelfing)
    initNoSelfing();

  if (_options._tree)
    initTree();

  initObligatoryLeaves();

  //if (!(_options._tree && _options._servin))
  initBounds();

  //if (_options._verbose)
  //{
  //  std::cout << "// Prob lower bound: " << log(_pData->getProbLowerBound())
  //            << "\t" << _pData->getProbLowerBound() << std::endl;
  //}
}

void IlpSolver::initPopExpr(IloExpr& expr) const
{
  expr.clear();
  for (size_t i = 0; i < _options._bound; i++)
  {
    for (size_t j = 0; j < _breakpoint.size(); j++)
    {
      expr += _lambda[i][j] * _N[j];
    }
  }
}

void IlpSolver::initBounds()
{
  //const size_t m = _pData->getNumberOfLoci();
  //const size_t n = _pData->getParents().size();
  //const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  LowerBound LB(_pData, false);

  // generations
  // fix generations of the backbone
  for (size_t i = 1; i < _options._fixedGen; i++)
  {
    _model.add(_r[_options._bound-i] == _options._fixedGen-i+1);
  }
  _model.add(_r[0] == 1);

  // use 1 <= gen < fixedGen for all non-root inner nodes
  for (size_t i = 0; i < _options._bound-1; i++)
  {
    if (!isBackboneNode(i))
    {
      _model.add(_r[i] >= 1);
      _model.add(_r[i] < _options._fixedGen);
    }
  }

  // bounds on pop
  IloExpr sumofpopsize(_env);
  initPopExpr(sumofpopsize);
  _model.add(sumofpopsize >= (int)LB.getPopLB());

  // at most the specified upper bond
  if (_options._upperBoundPop < std::numeric_limits<double>::max())
  {
    _model.add(sumofpopsize <= _options._upperBoundPop);
  }
  sumofpopsize.clear();
  sumofpopsize.end();

  // upper bound on objective
  if (_options._upperBoundObj < std::numeric_limits<double>::max())
  {
    _model.add(_obj <= _options._upperBoundObj);
  }
}

void IlpSolver::initSegments()
{
  _breakpoint.clear();
  double boundary = log(1-pow(1-_pData->getGamma(),1./_pData->getPopMax()));
  double seglength = (log(.5)-boundary)/_options._nof_segs;

  for (int i = 0; i <= _options._nof_segs; i++)
  {
    double val = boundary + i * seglength;
    _breakpoint.push_back(val);
  }
  _breakpoint.push_back(0);

  _N.clear();
  for (size_t j = 0; j < _breakpoint.size(); j++)
  {
    const double p = exp(_breakpoint[j]);
    const double pop = std::max(1., log(1-_pData->getGamma())/log(1-p));
    _N.push_back(pop);
  }
}

void IlpSolver::initParentalGenotypes()
{
  // initialize parental genotypes, store bits in _c
  const GenotypeSet& parents = _pData->getParents();
  const int m = _pData->getNumberOfLoci();

  int i=0;
  for (GenotypeSet::const_iterator it = parents.begin(); it != parents.end(); it++, i++)
  {
    const Genotype& parent = *it;

    _c[2*i] = IloBoolArray(_env, m);
    _c[2*i+1] = IloBoolArray(_env, m);
    for (int j=0; j<m; j++)
    {
      _c[2*i][j] = parent(m, 0, j);
      _c[2*i+1][j] = parent(m, 1, j);
    }
    _homozygous[i]=it->isHomozygous();
  }
}

void IlpSolver::initIdeotype(bool swapIdeotype)
{
  // initialize ideotype, store bits in _t1 and t2, assumed _t1=_t2
  const Genotype& ideotype = _pData->getIdeotype();
  const size_t m = _pData->getNumberOfLoci();

  // two targets - used only in multi-target case
  for (size_t j=0; j<m; j++)
  {
    _cs1[j] = ideotype(m, 0, j);
    _cs2[j] = ideotype(m, 1, j);
  }

  // fix target bits
  for (size_t j=0; j<m; j++)
    _model.add(_a[(2*_options._bound)-2][j]== (swapIdeotype ? _cs2[j] : _cs1[j]));

  if (!ideotype.isHomozygous())
  {
    for (size_t j=0; j<m; j++)
      _model.add(_a[(2*_options._bound)-1][j]== (swapIdeotype ? _cs1[j] : _cs2[j]));
  }
}

void IlpSolver::initChromosomesFromGenotypes()
{
  // initialize x variables and corresponding constraints
  const size_t n = _pData->getParents().size();
  const size_t g = _options._fixedGen;
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();

  // backbone chromosomes are the last 2*(g-1) chromosomes
  // + first two chromosomes

  // an x variable is defined as _x[k][i],
  // _x[k][i] = 1 iff chromosome k is obtained from node i
  // if k is backbone chromosome then 0 <= i < n+k/2
  // otherwise 0 <= i < n+k/2
  // OR _options._bound - _options._fixedGen <= i < _options._bound - 3

  for (size_t i = 0; i < _options._bound; i++)
  {
    // add offset if g is set and i is a non-backbone node
    int offset = getNrInnerPred(i);

    _x[2*i] = IloBoolVarArray(_env, n+offset);
    if (i < _options._bound - 1 || !homozygousIdeotype)
    {
      _x[2*i+1] = IloBoolVarArray(_env, n+offset);
    }

    for (size_t j = 0; j < n + offset; j++)
    {
      std::stringstream ss;
      ss << "x_" << 2*i << "_" << j;
      _x[2*i][j] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_x[2*i][j]);

      if (i < _options._bound - 1 || !homozygousIdeotype)
      {
        std::stringstream ss2;
        ss2 << "x_" << 2*i+1 << "_" << j;
        _x[2*i+1][j] = IloBoolVar(_env, ss2.str().c_str());
        _allVar.add(_x[2*i+1][j]);
      }
    }
  }

  // every chromosome has exactly one predecessor
  IloExpr sumofinedges1(_env);
  IloExpr sumofinedges2(_env);
  for (size_t i=0; i < _options._bound; i++)
  {
    int offset = getNrInnerPred(i);
    for (size_t j=0; j < n + offset; j++)
    {
      sumofinedges1 += _x[2*i][j];
      if (i<_options._bound-1 || !homozygousIdeotype)
        sumofinedges2 += _x[2*i+1][j];
    }
    _model.add(sumofinedges1 == 1);
    sumofinedges1.clear();

    if (i < _options._bound-1 || !homozygousIdeotype)
    {
      _model.add(sumofinedges2 == 1);
      sumofinedges2.clear();
    }
  }
  sumofinedges1.end();
  sumofinedges2.end();

  if (g > 1)
  {
    // out-degree for all *inner nodes* (except root node) must be at least 1
    IloExpr sumofoutedges(_env);
    for (size_t i = 0; i < _options._bound; i++)
    {
      // backbone nodes don't matter here, their out-degree is by definition >= 1
      if (!isBackboneNode(i))
      {
        // TODO: maybe upper bound for j can be tighter... DONE!
        // (REMEMBER: last node is only connected to the last but one)
        const size_t chromosomeUB = homozygousIdeotype ? 2 * _options._bound - 1 : 2 * _options._bound;
        for (size_t j = 2*(i+1); j < chromosomeUB; j++)
        {
          sumofoutedges += _x[j][n+i];
        }
        _model.add(sumofoutedges >= 1);
        sumofoutedges.clear();
      }
    }
    sumofoutedges.end();

    // fix the backbone
    // the backbone is a path of cardinality _options._fixedGen
    for (size_t h = 0; h < _options._fixedGen - 2; h++)
    {
      size_t i = _options._bound - h - 1;
      // we skip n+i-1, because we connect to n+i-1
      for (size_t j = 0; j < n+i-1; j++)
      {
        _model.add(_x[2*i][j] == 0);
      }
      _model.add(_x[2*i][n+i-1] == 1);
    }

    // connection from first backbone node to 2nd backbone node
    size_t i = _options._bound - _options._fixedGen + 1;
    for (size_t j = 0; j < n + i; j++)
    {
      if (j != n)
        _model.add(_x[2*i][j] == 0);
      else
        _model.add(_x[2*i][j] == 1);
    }
  }

  {
    // homozygous ideotype, last step is selfing:
    // REMOVE me, superfluous due to backbone constraints
    //for (int j = 0; j < n+_options._bound-2 ; j++)
    //{
    //  _model.add(_x[2*_options._bound-2][j] == 0);
    //}
    //_model.add(_x[2*_options._bound-2][n+_options._bound-2] == 1);
  }
}

void IlpSolver::initAllelesFromChromosomes()
{
  const size_t n = _pData->getParents().size();
  const size_t m = _pData->getNumberOfLoci();
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();

  // y-Vars, encodes bit origins
  // y[i][j] = { 0, if j-th bit of chromosome i originates 
  //                from the upper chromosome of the node defined by x
  //             1, if j-th bit of chromosome i originates
  //                from the lower chromosome of the node defined by x
  const size_t chromosomeUB = homozygousIdeotype ? 2 * _options._bound - 1 : 2 * _options._bound;
  for (size_t i = 0; i < chromosomeUB; i++)
  {
    _y[i] = IloBoolVarArray(_env, m);
    for (size_t j = 0; j < m; j++)
    {
      std::stringstream ss;
      ss << "y_" << i << "_" << j;
      _y[i][j] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_y[i][j]);
    }
  }

  // the bit variable
  for (size_t i = 0; i < chromosomeUB; i++)
  {
    _a[i] = IloBoolVarArray(_env, m);
    for (size_t t = 0; t < m; t++)
    {
      std::stringstream ss;
      ss << "a_" << i << "_" << t;
      _a[i][t] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_a[i][t]);
    }
  }

  // is locus heterozygous, 
  // (not defined for the target node, as it is homozygous by definition)
  const size_t genotypeUB = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  for (size_t i = 0; i < genotypeUB; i++)
  {
    _at[i] = IloBoolVarArray(_env, m);
    for (size_t t = 0; t < m; t++)
    {
      std::stringstream ss;
      ss << "at_" << i << "_" << t;
      _at[i][t] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_at[i][t]);
    }
  }

  for (size_t i = 0; i < genotypeUB; i++)
  {
    for (size_t j = 0; j < m; j++)
    {
      // TODO reicht das? keine <= Ungleichung notwendig?
      // Ja das reicht! Being homozygous is beneficial
      _model.add(_at[i][j] >= _a[2*i][j] - _a[2*i+1][j]);
      _model.add(_at[i][j] >= _a[2*i+1][j] - _a[2*i][j]);
    }
  }

  // is node heterozygous
  for (size_t i = 0; i < genotypeUB; i++)
  {
    std::stringstream ss;
    ss << "h_" << i;
    _h[i] = IloBoolVar(_env, ss.str().c_str());
    _allVar.add(_h[i]);
    for (size_t j = 0; j < m; j++)
      _model.add(_h[i] >= _at[i][j]);
  }


  // selfing lemma:
  // if a node i is homozygous (i.e. _het[i]==0), 
  // then it's two predecessors must be the same
  for (size_t i = 0; i < genotypeUB; i++)
  {
    int offset = getNrInnerPred(i);
    for (size_t j = 0; j < n + offset; j++)
    {
      _model.add(_x[2*i][j] - _x[2*i+1][j] <= _h[i]);
      _model.add(_x[2*i+1][j] - _x[2*i][j] <= _h[i]);
    }
  }

  // _g[i][j][t] = 1 iff 1-allele at locus j from chr i comes from chr t
  for (size_t i = 0; i < chromosomeUB; i++)
  {
    _g[i] = BoolVarMatrix(_env, m);

    size_t offset = 2*getNrInnerPred(i/2);
    for (size_t j = 0; j < m; j++)
    {
      _g[i][j] = IloBoolVarArray(_env, 2*n + offset);

      for (size_t t = 0; t < 2*n + offset; t++)
      {
        std::stringstream ss;
        ss << "g_";
        ss << i << "_" << j << "_" << t;
        _g[i][j][t] = IloBoolVar(_env, ss.str().c_str());
        _allVar.add(_g[i][j][t]);
      }
    }
  }

  IloExpr sumof1dep(_env);
  for (size_t i = 0; i < chromosomeUB; i++)
  {
    for (size_t j = 0; j < m; j++)
    {
      // inner nodes
      size_t offset = 2*getNrInnerPred(i/2);
      for (size_t l = 0; l < offset; l++)
      {
        _model.add(_g[i][j][2*n+l] <= _x[i][n + l/2]);
        _model.add(_g[i][j][2*n+l] <= _a[getAbsChromosome(i/2, l)][j]);
        if ((l%2) == 0)
        {
          _model.add(_g[i][j][2*n+l] <= 1-_y[i][j]);
          _model.add(_g[i][j][2*n+l] >= _x[i][l/2+n] + _a[getAbsChromosome(i/2,l)][j] + (1-_y[i][j])-2);
        }
        else
        {
          _model.add(_g[i][j][2*n+l] <= _y[i][j]);
          _model.add(_g[i][j][2*n+l] >= _x[i][l/2+n] + _a[getAbsChromosome(i/2,l)][j] + _y[i][j]-2);
        }
      }

      // leaves
      for (size_t l = 0; l < 2*n; l++)
      {
        _model.add(_g[i][j][l] <= _c[l][j] * _x[i][l/2]);
        if ((l%2) == 0)
        {
          _model.add(_g[i][j][l] <= _c[l][j] * (1-_y[i][j]));
          _model.add(_g[i][j][l] >= _c[l][j] * (_x[i][l/2]+(1-_y[i][j])-1));
        }
        else
        {
          _model.add(_g[i][j][l] <= _c[l][j] * _y[i][j]);
          _model.add(_g[i][j][l] >= _c[l][j] * (_x[i][l/2]+_y[i][j]-1));
        }
      }  
      for (size_t l=0; l < 2*n + offset; l++)
      {
        sumof1dep += _g[i][j][l];
      }
      _model.add(sumof1dep == _a[i][j]);
      sumof1dep.clear();
    }
  }
  sumof1dep.end();
}

void IlpSolver::initUsefulCross()
{ 
  const size_t m = _pData->getNumberOfLoci();
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const size_t chromosomeUB = homozygousIdeotype ? 2 * _options._bound - 1 : 2 * _options._bound;

  IloExpr sumofcollected(_env);
  for (size_t i=0; i < chromosomeUB; i++)
  {
    for (size_t j=0; j<m-2; j++)
    {
      for (size_t l=j+1; l<m-1; l++)
      {
        for (size_t t=j+1; t<l+1; t++)
          sumofcollected+=_a[i][t];
        _model.add(_d[i][j]+_d[i][l]-sumofcollected<=1);
        sumofcollected.clear();
      }
    }

    for (size_t j=0; j<m-1; j++)
    {
      for (size_t t=0; t<j+1; t++)
        sumofcollected+=_a[i][t];
      _model.add(_d[i][j]-sumofcollected<=0);
      sumofcollected.clear();
    }      
    for (size_t j=0; j<m-1; j++)
    {
      for (size_t t=j+1; t<m; t++)
        sumofcollected+=_a[i][t];
      _model.add(_d[i][j]-sumofcollected<=0);
      sumofcollected.clear();
    }
  }
  sumofcollected.end();
}

void IlpSolver::initNoSelfing()
{
  const size_t n = _pData->getParents().size();
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const size_t genotypeUB = homozygousIdeotype ? _options._bound - 1 : _options._bound;

  for (size_t i = 0; i < genotypeUB; i++)
  {
    if (homozygousIdeotype && i == genotypeUB - 1) continue;

    int offset = getNrInnerPred(i);
    for (size_t j = 0; j < n + offset; j++)
      _model.add(_x[2*i][j]+_x[2*i+1][j]<=1);
  }
}

void IlpSolver::initTree()
{
  abort();
  /*
  // TODO: fix this!
  const int m = _pData->getNumberOfLoci();
  const int n = _pData->getParents().size();
  const int g = _options._fixedGen;

  // every node is used exactly once as parent (=tree, out-degree = 1)
  // in other words: every node has only one successor
  IloExpr sumofoutedges(_env);

  // leaves are used at most once
  for (int i=0; i<n; i++)
  {
    for (int j=0; j < (2 * _options._bound)-1; j++)
    {
      sumofoutedges+=_x[j][i];
    }
    // TODO maybe at most 1 here? DONE!
    _model.add(sumofoutedges<=1);
    sumofoutedges.clear();
  }

  if (g > 1)
  {
    for (int i=0; i<_options._bound-g; i++)
    {
      // successors of a non-backbone inner node i 
      // are the nodes coming after i in the topological ordering
      // except last two nodes
      for (int j=2*(i+1); j < (2 * _options._bound)-1; j++)
      {
        sumofoutedges+=_x[j][n+i];
      }
      _model.add(sumofoutedges==1).setName("outdeg");
      sumofoutedges.clear();
    }
    for (int i=_options._bound-g; i<_options._bound-2; i++)
    {
      // the successor of a backbone inner node i is the immediate successor (i+1)
      // everything else cannot be a successor

      // backbone nodes only have non-backbone inner nodes as successor as well 
      // as the next backbone node
      for (int j=0; j < 2*(_options._bound-g); j++)
      {
        sumofoutedges+=_x[j][n+getRelNode(j/2, i)];
      }
      // TODO what about the succeeding backbone nodes
      _model.add(sumofoutedges==0).setName("outdegbb");
      sumofoutedges.clear();
    }
  }
  else
  {
    for (int i=0; i<_options._bound-1; i++)
    {
      for (int j=2*(i+1); j < (2 * _options._bound)-1; j++)
      {
        sumofoutedges+=_x[j][n+i];
      }
      _model.add(sumofoutedges==1);
      sumofoutedges.clear();
    }
  }
  sumofoutedges.end();

  // uniquely covered stuff (t1 assumed to be equal to t2)
  // cover[p] how many leaves cover the target allele at p
  std::vector<int> cover(m, 0);
  for (int p = 0; p < m; p++)
  {
    for (int i = 0; i < n; i++)
    {
      if (_cs1[p] == _c[2*i][p] || _cs1[p] == _c[2*i+1][p])
      {
        cover[p]++;
      }
    }
  }

  // _uc[i][p] = 1 iff locus p at node i contains a target allele that cannot be lost
  // _ucx[k][p][i] chromosome k depends (x) on node i that has a unique 1 in column p (uc)
  // i.e.: _ucx[k][p][i] = _x[k][i]*_uc[i][p]
  for (int i = 0; i < _options._bound; i++)
  {
    _uc[i] = IloBoolVarArray(_env, m);

    _ucx[2*i] = BoolVarMatrix(_env, m);
    if (i != _options._bound - 1)
      _ucx[2*i+1] = BoolVarMatrix(_env, m);

    for (int p = 0; p < m; p++)
    {
      _uc[i][p] = IloBoolVar(_env);
      _allVar.add(_uc[i][p]);

      int offset = g <= 1 || i >= _options._bound-g ? 0 : g-2;
      
      _ucx[2*i][p] = IloBoolVarArray(_env, n+i+offset);
      for (int j = n; j < n+i+offset; j++)
      {
        _ucx[2*i][p][j-n] = IloBoolVar(_env);
        _allVar.add(_ucx[2*i][p][j-n]);
      }
      
      if (i != _options._bound - 1)
      {
        _ucx[2*i+1][p] = IloBoolVarArray(_env, n+i+offset);
        for (int j = n; j < n+i+offset; j++)
        {
          _ucx[2*i+1][p][j-n] = IloBoolVar(_env);
          _allVar.add(_ucx[2*i+1][p][j-n]);
        }
      }
    }
  }

  // _ucx[k][p][i] = _x[k][i]*_uc[i][p]
  // we are going to model:
  // uc[j][p] = sum_{i in pred(j)} (ucx[2j][p][i] + ucx[2j+1][p][i])
  IloExpr sum1(_env);
  IloExpr sum2(_env);
  for (int i = 0; i < _options._bound; i++)
  {
    for (int p = 0; p < m; p++)
    {
      // first the leaves
      for (int j = 0; j < n; j++)
      {
        if (cover[p] == 1 && _c[2*j][p] == _cs1[p])
        {
          sum1+=_x[2*i][j];
        }
        if (i < _options._bound-1 && cover[p] == 1 && _c[2*j+1][p] == _cs1[p])
        {
          sum2+=_x[2*i+1][j];
        }
      }
      // then the internal nodes
      int offset = g <= 1 || i >= _options._bound-g ? 0 : g-2;
      for (int j = n; j < n+i+offset; j++)
      {
        sum1+=_ucx[2*i][p][j-n];

        _model.add(_ucx[2*i][p][j-n] <= _uc[getAbsNode(i,j-n)][p]);
        _model.add(_ucx[2*i][p][j-n] <= _x[2*i][j]);
        _model.add(_ucx[2*i][p][j-n] >= _uc[getAbsNode(i,j-n)][p] + _x[2*i][j] - 1);

        if (i < _options._bound-1)
        {
          sum2+=_ucx[2*i+1][p][j-n];

          _model.add(_ucx[2*i+1][p][j-n] <= _uc[getAbsNode(i,j-n)][p]);
          _model.add(_ucx[2*i+1][p][j-n] <= _x[2*i+1][j]);
          _model.add(_ucx[2*i+1][p][j-n] >= _uc[getAbsNode(i,j-n)][p] + _x[2*i+1][j] - 1);
        }
      }
      _model.add(_uc[i][p] == sum1 + sum2);

      if (_options._servin)
      {
        _model.add(_a[2*i][p] == sum1);
        if (i < _options._bound-1)
        {
          _model.add(_a[2*i+1][p] == sum2);
          // homozygous => 0,0
          _model.add(_a[2*i][p] + _a[2*i+1][p] <= _at[i][p]);
        }
      }
      else
      {
        // if locus p at chromosome 2i is unique, we must keep it
        _model.add(_a[2*i][p] >= sum1);
        if (i < _options._bound-1)
          // if locus p at chromosome 2i+1 is unique, we must keep it
          _model.add(_a[2*i+1][p] >= sum2);
      }
      sum1.clear();
      sum2.clear();
    }
  }
  sum1.end();
  sum2.end();*/
}

void IlpSolver::initObligatoryLeaves()
{
  const size_t m = _pData->getNumberOfLoci();
  const size_t n = _pData->getParents().size();
  const GenotypeSet& parents = _pData->getParents();
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();

  // outdegree of obligatory leaves is at least one
  bool* obligatory = new bool[n];
  for (size_t i=0; i<n; i++)
    obligatory[i]=false;
 
  for (size_t j=0; j<m; j++)
  {
    size_t i=0;
    bool unique_covered = false;
    size_t obl_parent = 0;
    for (GenotypeSet::const_iterator it = parents.begin(); it != parents.end(); it++)
    {
      if (_c[2*i][j] || _c[2*i+1][j])
      {
        if (unique_covered) {
          unique_covered=false;
          break;
        }
        else {
          unique_covered=true;
          obl_parent = i;
        }
      }
      i++;
    }
    if (unique_covered)
    {     
      obligatory[obl_parent]=true;
    }
  }

  const size_t chromosomeUB = homozygousIdeotype ? 2 * _options._bound - 1 : 2 * _options._bound;
  IloExpr sumofoutedges(_env);
  for (size_t i=0; i<n; i++)
  {
    if (!obligatory[i])
      continue;
    for (size_t j=0; j < chromosomeUB; j++)
    {    
      sumofoutedges+=_x[j][i];
    }
    _model.add(sumofoutedges>=1);       
    sumofoutedges.clear();
  }    
  sumofoutedges.end();
  delete obligatory;
}

void IlpSolver::initGen()
{
  const size_t n = _pData->getParents().size();
  const size_t bound = _options._bound;
  const int g = static_cast<int>(_options._fixedGen);

  for (size_t i = 0; i < bound; i++)
  {
    std::stringstream ss;
    ss << "r_" << i;
    _r[i].setName(ss.str().c_str());
    _allVar.add(_r[i]);

    size_t offset = getNrInnerPred(i);
    for (size_t j = 0; j < offset; j++)
    {
      _model.add(_r[i] >= _r[getAbsNode(i,j)] + 1 - (g*(1-_x[2*i][n+j])));
      if (i < bound-1)
        _model.add(_r[i] >= _r[getAbsNode(i,j)] + 1 - (g*(1-_x[2*i+1][n+j])));
    }
  }
}

void IlpSolver::initPop()
{
  const size_t n = _pData->getParents().size();
  const size_t m = _pData->getNumberOfLoci();
  const DoubleMatrix& RM = _pData->getRM();
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();

  // product var for heterozygous and x
  // _hx[k][i] = _hx[i]*_x[k][i]
  for (size_t i = 0; i < _options._bound; i++)
  {
    int offset = getNrInnerPred(i);
    _hx[2*i] = IloBoolVarArray(_env, n + offset);
    if (i <_options._bound - 1 || !homozygousIdeotype)
    {
      _hx[2*i+1] = IloBoolVarArray(_env, n + offset);
    }

    // leaves
    for (size_t j = 0; j < n; j++)
    { 
      std::stringstream ss;
      ss << "hx_" << 2*i << "_" << j;
      _hx[2*i][j] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_hx[2*i][j]);

      if (i < _options._bound - 1 || !homozygousIdeotype)
      {
        std::stringstream ss2;
        ss2 << "hx_" << 2*i + 1 << "_" << j;
        _hx[2*i+1][j] = IloBoolVar(_env, ss2.str().c_str());
        _allVar.add(_hx[2*i+1][j]);
      }
      if (_homozygous[j])
      {
        _model.add(_hx[2*i][j]==0);
        if (i<_options._bound-1 || !homozygousIdeotype)
        {
          _model.add(_hx[2*i+1][j]==0);
        }
      }
      else
      {
        _hx[2*i][j] = _x[2*i][j];
        if (i<_options._bound-1 || !homozygousIdeotype)
        {
          _hx[2*i+1][j] = _x[2*i+1][j];
        }
      }
    }

    // inner nodes
    for (size_t j = n; j < n + offset; j++)
    {
      std::stringstream ss;
      ss << "hx_" << 2*i << "_" << j;
      _hx[2*i][j] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_hx[2*i][j]);

      if (i < _options._bound - 1 || !homozygousIdeotype)
      {
        std::stringstream ss2;
        ss2 << "hx_" << 2*i + 1 << "_" << j;
        _hx[2*i+1][j] = IloBoolVar(_env, ss2.str().c_str());
        _allVar.add(_hx[2*i+1][j]);
      }

      _model.add(_hx[2*i][j]<=_h[getAbsNode(i,j-n)]);
      _model.add(_hx[2*i][j]<=_x[2*i][j]);
      _model.add(_hx[2*i][j]>=_h[getAbsNode(i,j-n)]+_x[2*i][j]-1);

      if (i < _options._bound - 1 || !homozygousIdeotype)
      {
        _model.add(_hx[2*i+1][j]<=_h[getAbsNode(i,j-n)]);
        _model.add(_hx[2*i+1][j]<=_x[2*i+1][j]);
        _model.add(_hx[2*i+1][j]>=_h[getAbsNode(i,j-n)]+_x[2*i+1][j]-1);
      }
    }
  }

  // did jth bit in string i result from a crossover event?
  const size_t chromosomeUB = homozygousIdeotype ? 2 * _options._bound - 1 : 2 * _options._bound;
  for (size_t i = 0; i < chromosomeUB; i++)
  {
    _d[i] = IloBoolVarArray(_env, m-1);
    for (size_t j = 0; j < m - 1; j++)
    {
      std::stringstream ss;
      ss << "d_" << i << "_" << j;
      _d[i][j] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_d[i][j]);
    }
  }

  for (size_t i = 0; i < chromosomeUB; i++)
  {
    for (size_t j = 1; j < m; j++)
    {
      _model.add(_d[i][j-1] >= _y[i][j]-_y[i][j-1]);
      _model.add(_d[i][j-1] >= _y[i][j-1]-_y[i][j]);
    }
  }

  // b_{i,p,q} = 1 iff loci p and q of node i are heterozygous
  // and everything in between is homozygous
  const size_t genotypeUB = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  IloExpr sumofr(_env);
  for (size_t i = 0; i < genotypeUB; i++)
  {
    _b[i] = BoolVarMatrix(_env, m - 1);
    for (size_t p = 0; p < m - 1; p++)
    {
      _b[i][p] = IloBoolVarArray(_env, m-p-1);
      for (size_t q = p + 1; q < m; q++)
      {
        std::stringstream ss;
        ss << "b_" << i << "_" << p << "_" << q-p-1;
        _b[i][p][q-p-1] = IloBoolVar(_env, ss.str().c_str());
        _allVar.add(_b[i][p][q-p-1]);
        _model.add(_b[i][p][q-p-1]<=_at[i][p]);
        _model.add(_b[i][p][q-p-1]<=_at[i][q]);
        for (size_t r = p + 1; r < q; r++)
        {
          _model.add(_b[i][p][q-p-1]<=1-_at[i][r]);
          sumofr += 1-_at[i][r];
        }
        _model.add(_b[i][p][q-p-1] >= _at[i][p] + _at[i][q] + sumofr - (q-p));
        sumofr.clear();
      }
    }
  }
  sumofr.end();

  // product var of _b and _x:
  // _bx[k][i][p][q] = _b[i][p][q]*_x[k][i]
  for (size_t k = 0; k < chromosomeUB; k++)
  {
    int offset = getNrInnerPred(k/2);
    _bx[k] = BoolVar3Matrix(_env, n + offset);
    
    // leaves
    for (size_t i = 0; i < n; i++)
    {
      _bx[k][i] = BoolVarMatrix(_env, m-1);
      for (size_t p = 0; p < m-1; p++)
      {
        _bx[k][i][p] = IloBoolVarArray(_env, m-p-1);
        for (size_t q = p+1; q < m; q++)
        {
          std::stringstream ss;
          ss << "bx_" << k << "_" << i << "_" << p << "_" << q-p-1;
          _bx[k][i][p][q-p-1] = IloBoolVar(_env, ss.str().c_str());
          _allVar.add(_bx[k][i][p][q-p-1]);
          if (_c[2*i][p] == _c[2*i+1][p] ||
              _c[2*i][q] == _c[2*i+1][q] ||
              !isHomozygousBlock(i,p,q))
          {
            // _h[i][p][q] == 0
            _model.add(_bx[k][i][p][q-p-1] == 0);
          }
          else
          {
            _model.add(_bx[k][i][p][q-p-1] == _x[k][i]);
          }
        }
      }
    }

    // inner nodes
    for (size_t i = n; i < n + offset; i++)
    {
      _bx[k][i] = BoolVarMatrix(_env, m-1);
      for (size_t p = 0; p < m-1; p++)
      {
        _bx[k][i][p] = IloBoolVarArray(_env, m-p-1);
        for (size_t q = p+1; q < m; q++)
        {
          std::stringstream ss;
          ss << "bx_" << k << "_" << i << "_" << p << "_" << q-p-1;
          _bx[k][i][p][q-p-1] = IloBoolVar(_env, ss.str().c_str());
          _allVar.add(_bx[k][i][p][q-p-1]);
          _model.add(_bx[k][i][p][q-p-1] <= _b[getAbsNode(k/2,i-n)][p][q-p-1]);
          _model.add(_bx[k][i][p][q-p-1] <= _x[k][i]);
          _model.add(_bx[k][i][p][q-p-1] >= _b[getAbsNode(k/2,i-n)][p][q-p-1] + _x[k][i] - 1);
        }
      }
    }
  }

  // var z: z[k][p][q] = 1 iff the alleles at loci p and q of chr k
  // are the result of a crossover
  IloExpr sumofr2(_env);
  for (size_t k = 0; k < chromosomeUB; k++)
  {
    _z[k] = IntVarMatrix(_env, m-1);

    for (size_t p = 0; p < m-1; p++)
    {
      _z[k][p] = IloIntVarArray(_env, m-p-1);
      for (size_t q = p+1; q < m; q++)
      {
        std::stringstream ss;
        ss << "z_" << k << "_" << p << "_" << q-p-1;
        _z[k][p][q-p-1] = IloIntVar(_env, 0, m-1, ss.str().c_str());
        _allVar.add(_z[k][p][q-p-1]);
        for (size_t r = p + 1; r <= q; r++)
        {
          sumofr2 += _d[k][r-1]; // mind the r-1
        }
        _model.add(_z[k][p][q-p-1] == sumofr2);
        sumofr2.clear();
      }
    }
  }
  sumofr2.end();

  // product var of _b, _x and _z:
  // _bxz[k][i][p][q] = _b[i][p][q]*_x[k][i]*_z[k][p][q]
  for (size_t k=0; k < chromosomeUB; k++)
  {
    int offset = getNrInnerPred(k/2);
    _bxz[k] = IntVar3Matrix(_env, n + offset);
    
    for (size_t i = 0; i < n; i++)
    {
      _bxz[k][i] = IntVarMatrix(_env, m - 1);
      for (size_t p = 0; p < m - 1; p++)
      {
        _bxz[k][i][p] = IloIntVarArray(_env, m-p-1);
        for (size_t q = p+1; q < m; q++)
        {
          std::stringstream ss;
          ss << "bxz_" << k << "_" << i << "_" << p << "_" << q-p-1;
          _bxz[k][i][p][q-p-1] = IloIntVar(_env, 0, m-1, ss.str().c_str());
          _allVar.add(_bxz[k][i][p][q-p-1]);
          if (_c[2*i][p] == _c[2*i+1][p] ||
              _c[2*i][q] == _c[2*i+1][q] ||
              !isHomozygousBlock(i,p,q))
          {
            // _b[i][p][q] == 0
            _model.add(_bxz[k][i][p][q-p-1] == 0);
          }
          else
          {
            _model.add(_bxz[k][i][p][q-p-1] <= _x[k][i] * static_cast<int>(m-1));
            _model.add(_bxz[k][i][p][q-p-1] <= _z[k][p][q-p-1]);
            _model.add(_bxz[k][i][p][q-p-1] >= _x[k][i] + (1./(m-1)) * _z[k][p][q-p-1] - 1);
          }
        }
      }
    }
    for (size_t i = n; i < n + offset; i++)
    {
      _bxz[k][i] = IntVarMatrix(_env, m-1);
      for (size_t p = 0; p < m-1; p++)
      {
        _bxz[k][i][p] = IloIntVarArray(_env, m-p-1);
        for (size_t q = p+1; q < m; q++)
        {
          std::stringstream ss;
          ss << "bxz_" << k << "_" << i << "_" << p << "_" << q-p-1;
          _bxz[k][i][p][q-p-1] = IloIntVar(_env, 0, m-1, ss.str().c_str());
          _allVar.add(_bxz[k][i][p][q-p-1]);
          // TODO, maybe better in terms of b and x
          _model.add(_bxz[k][i][p][q-p-1] <= _bx[k][i][p][q-p-1] * static_cast<int>(m-1));
          _model.add(_bxz[k][i][p][q-p-1] <= _z[k][p][q-p-1]);
          _model.add(_bxz[k][i][p][q-p-1] >= _bx[k][i][p][q-p-1] + (1./(m-1)) * _z[k][p][q-p-1] - 1);
        }
      }
    }
  }

  // separable p var
  IloExpr sumoflogprobs(_env);
  for (size_t j = 0; j < _options._bound; j++)
  {
    std::stringstream ss;
    ss << "p_" << j;
    _p[j].setName(ss.str().c_str());
    _allVar.add(_p[j]);

    size_t offset = getNrInnerPred(j);
    for (size_t i = 0; i < n + offset; i++)
    {
      if (i < n && !_homozygous[i])
      {
        // factor of 0.5 if parent is heterozygous
        sumoflogprobs += log(0.5) * _x[2*j][i];
        if (j < _options._bound - 1 || !homozygousIdeotype)
          sumoflogprobs += log(0.5)*_x[2*j+1][i];
      }
      else if (i >= n)
      {
        // factor of 0.5 of parent is heterozygous
        sumoflogprobs += log(0.5)*_hx[2*j][i];
        if (j < _options._bound - 1 || !homozygousIdeotype)
          sumoflogprobs+=log(0.5)*_hx[2*j+1][i];
      }

      for (size_t p = 0; p < m - 1; p++)
      {
        for (size_t q = p + 1; q < m; q++)
        {
          double r_pq = RM[p][q];
          sumoflogprobs += log(1-r_pq) * _bx[2*j][i][p][q-p-1];
          sumoflogprobs += log(r_pq/(1.-r_pq)) * _bxz[2*j][i][p][q-p-1];

          if (j < _options._bound - 1 || !homozygousIdeotype)
          {
            sumoflogprobs += log(1-r_pq) * _bx[2*j+1][i][p][q-p-1];
            sumoflogprobs += log(r_pq/(1.-r_pq)) * _bxz[2*j+1][i][p][q-p-1];
          }
        }
      }
    }
    if (j < _options._bound - 1 || !homozygousIdeotype)
    {
      _model.add(_p[j] == sumoflogprobs);
      if (_options._boundNmax)
        _model.add(sumoflogprobs >= log(_pData->getProbLowerBound()));
    }
    else if (homozygousIdeotype)
    {
      // last crossing is a selfing (assumption: homozygous ideotype)
      _model.add(_p[j] == sumoflogprobs + sumoflogprobs);
      if (_options._boundNmax)
        _model.add(sumoflogprobs + sumoflogprobs >= log(_pData->getProbLowerBound()));
    }
    sumoflogprobs.clear();
  }
  sumoflogprobs.end();
}

void IlpSolver::initObj()
{
  IloExpr sumoflambdas(_env);
  IloExpr sumoflambdas1(_env);
  for (size_t i = 0; i < _options._bound; i++)
  {
    _lambda[i] = IloNumVarArray(_env, _breakpoint.size(), 0., 1., IloNumVar::Float);

    for (size_t j = 0; j < _breakpoint.size(); j++)
    {
      std::stringstream ss;
      ss << "lambda_" << i << "_" << j;
      if (j == _breakpoint.size() - 1)
      {
        _lambda[i][j] = IloBoolVar(_env);
      }
      _lambda[i][j].setName(ss.str().c_str());
      _allVar.add(_lambda[i][j]);

      sumoflambdas1 += _breakpoint[j] * _lambda[i][j];
      sumoflambdas += _lambda[i][j];
    }

    _model.add(sumoflambdas1 == _p[i]);
    _model.add(sumoflambdas == 1);

    sumoflambdas.clear();
    sumoflambdas1.clear();
  }
  sumoflambdas.end();
  sumoflambdas1.end();

  _obj.clear();
  initPopExpr(_obj);
}

void IlpSolver::constructDAG()
{
  const GenotypeSet& parents = _pData->getParents();
  const int n = parents.size();

  GenotypePointerVector genotypeVector;
  std::vector<Node> constructed(n + _options._bound, lemon::INVALID);

  // reconstruct parental genotypes
  for (GenotypeSet::const_iterator it = parents.begin();
    it != parents.end(); it++)
  {
    genotypeVector.push_back(new Genotype(*it));
  }

  int targetIdx = n + _options._bound - 1;

  // reconstruct genotypes of inner nodes
  for (size_t i = 0; i < _options._bound-1; i++)
  {
    genotypeVector.push_back(new Genotype(parseGenotype(i)));
  }
  genotypeVector.push_back(new Genotype(_pData->getIdeotype()));

  // reconstruct arcs and nodes
  bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  for (size_t i = 0; i < _options._bound; i++)
  {
    if (!genotypeVector[n+i])
      continue;

    int p1 = -1, p2 = -1;
    size_t offset = getNrInnerPred(i);
    for (size_t j = 0; j < n + offset; j++)
    {
      double x_val = _pCplex->getValue(_x[2*i][j]);
      if (_tol.nonZero(x_val))
      {
        assert(p1 == -1);
        if (j>=n+i)
          p1 = n+getAbsNode(i, j-n);
        else
          p1 = j;
      }
      if (homozygousIdeotype)
      {
        if (i != _options._bound - 1 && _tol.nonZero(_pCplex->getValue(_x[2*i+1][j])))
        {
          assert(p2 == -1);
          if (j>=n+i)
            p2 = n+getAbsNode(i, j-n);
          else
            p2 = j;
        }
        else if (i == _options._bound - 1)
          p2 = p1;
      }
      else
      {
        double x_val2 = _pCplex->getValue(_x[2*i+1][j]);
        if (_tol.nonZero(x_val2))
        {
          assert(p2 == -1);
          if (j>=n+i)
            p2 = n+getAbsNode(i, j-n);
          else
            p2 = j;
        }
      }
    }

    assert(p1 != -1 || p2 != -1);
    Node node = constructed[n+i];
    if (node == lemon::INVALID) node = _G.addNode();
    _genotype[node] = *genotypeVector[n+i];
    _pop[node] = parsePop(i);
    double prob1 = parseProb1(i);
    double prob2 = parseProb2(i);
    _prob[node] = prob1 + prob2;
    _prob1[node] = prob1;
    _prob2[node] = prob2;
    _gen[node] = parseGen(i);
    _idx[node] = i;
    constructed[n+i] = node;
    
    if (p1 != -1)
    {
      Node nodeP1 = constructed[p1];
      if (nodeP1 == lemon::INVALID)
      {
        nodeP1 = _G.addNode();
        _genotype[nodeP1] = *genotypeVector[p1];
        _gen[nodeP1] = 0;
        _pop[nodeP1] = 0;
        _prob[nodeP1] = 0;
        _idx[nodeP1] = -1;
        constructed[p1] = nodeP1;
      }

      // add new arc
      Arc a1 = _G.addArc(nodeP1, node);
      _upperChr[a1] = true;
    }

    if (p2 != -1)
    {
      Node nodeP2 = constructed[p2];
      if (nodeP2 == lemon::INVALID)
      {
        nodeP2 = _G.addNode();
        _genotype[nodeP2] = *genotypeVector[p2];
        _gen[nodeP2] = 0;
        _pop[nodeP2] = 0;
        _prob[nodeP2] = 0;
        _idx[nodeP2] = -1;
        constructed[p2] = nodeP2;
      }

      // add new arc
      Arc a2 = _G.addArc(nodeP2, node);
      _upperChr[a2] = false;
    }
  }

  _targetNode = constructed[targetIdx];

  BoolNodeMap visited(_G, false);
  updateAttributesDAG(_targetNode, visited, false);

  // clean-up
  for (GenotypePointerVector::iterator it = genotypeVector.begin(); 
    it != genotypeVector.end(); it++)
  {
    delete *it;
  }
}
