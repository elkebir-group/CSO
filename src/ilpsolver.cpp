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
  , _idx(_G)
  , _env()
  , _model(_env)
  , _allVar(_env)
  , _obj(_env)
  , _pCplex(NULL)
  , _breakpoint()
  , _t1(_env, pData->getNumberOfLoci())
  , _t2(_env, pData->getNumberOfLoci())
  , _c(_env, 2 * pData->getParents().size())
  , _homozygous(_env, pData->getParents().size())
  , _x(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _xx(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _input_bit(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _y(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _p(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _bt(_env, pData->isIdeotypeHomozygous() ? options._bound-1 : options._bound)
  , _het(_env, pData->isIdeotypeHomozygous() ? options._bound-1 : options._bound)
  , _hetx(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _d(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _zb(_env, options._bound, -IloInfinity, 0.0, IloNumVar::Float)
  , _lambda(_env, options._bound)
  , _r(_env, options._bound, 1, options._fixedGen)
  , _h(_env, pData->isIdeotypeHomozygous() ? options._bound-1 : options._bound)
  , _hx(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _hxe(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _e(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
  , _uc(_env, options._bound)
  , _ucx(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-1 : 2 * options._bound)
{
}

IlpSolver::~IlpSolver()
{
  delete _pCplex;
  _obj.end();
  _t1.end();
  _t2.end();
  _c.end();
  _homozygous.end();
  _x.end();
  _xx.end();
  _input_bit.end();
  _y.end();
  _p.end();
  _bt.end();
  _het.end();
  _hetx.end();
  _d.end();
  _zb.end();
  _lambda.end();
  _r.end();
  _h.end();
  _hx.end();
  _hxe.end();
  _e.end();
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

  //char buf[1024];
  //sprintf(buf, "letssee-c%d-g%d.lp", _options._bound, _options._fixedGen);
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
    std::cout << "// Objective value: " << _pCplex->getObjValue() << std::endl;
    //printInnerNodes();
    //printH();
    std::cout << "//" << std::endl;
    printX();
    std::cout << "//" << std::endl;
    printXX();
    //std::cout << "//" << std::endl;
    //printE();
    //std::cout << "//" << std::endl;
    //printHX();
    //std::cout << "//" << std::endl;
    //printHXE();
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

  initObjectiveGen();

  // TODO do we need "|| _pData->getCostCrossover()" ??
  if (_pData->getCostNode() || _pData->getCostCrossover())
    initObjectiveCross();

  if (_pData->getCostCrossover())
    initObjectivePop();

  _model.add(IloMinimize(_env, _obj));

  if (_options._usefulCross)
    initUsefulCross();

  if (_options._diffParents)
    initDiffParents();

  if (_options._tree)
    initTree();

  initObligatoryLeaves();

  //if (!(_options._tree && _options._servin))
  initBounds();
}

void IlpSolver::initBounds()
{
  const int m = _pData->getNumberOfLoci();
  const int n = _pData->getParents().size();
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  LowerBound LB(_pData, false);

  // generations
  // fix generations of the backbone
  for (int i = 1; i < _options._fixedGen; i++)
  {
    _model.add(_r[_options._bound-i] == _options._fixedGen-i+1);
  }
  _model.add(_r[0] == 1);

  // use 1 <= gen < fixedGen for all non-root inner nodes
  for (int i = 0; i < _options._bound-1; i++)
  {
    if (!isBackboneNode(i))
    {
      _model.add(_r[i] >= 1);
      _model.add(_r[i] < _options._fixedGen);
    }
  }

  // bounds on pop
  IloExpr sumofpopsize(_env);
  for (int i=0; i < _options._bound; i++)
  {
    for (int j=0; j < _breakpoint.size(); j++)
    {
      if (j == _breakpoint.size() - 1)
      {
        assert(_breakpoint[j] == 0 && _pData->getGamma() >= .5);
        sumofpopsize+=_lambda[i][j];
      }
      else
      {
        sumofpopsize+=(log(1-_pData->getGamma())/log(1-exp(_breakpoint[j])))*_lambda[i][j];
      }
    }
  }
  // at least the lower bound, TODO: give lemma here
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
  double l_boundary = log(1-pow(1-_pData->getGamma(),1./_pData->getPopMax()));
  double seglength = -1*(l_boundary-log(.5))/_options._nof_segs;

  for (int i=0; i<_options._nof_segs; i++)
    _breakpoint.push_back(l_boundary+(i*seglength));    

  // add default breakpoints
  _breakpoint.push_back(log(.5));
  _breakpoint.push_back(0);
}

void IlpSolver::initParentalGenotypes()
{
  // initialize parental genotypes, store bits in _c
  const GenotypeSet& parents = _pData->getParents();
  const int m = _pData->getNumberOfLoci();

  int i=0;
  for (GenotypeSet::const_iterator it = parents.begin(); it != parents.end(); it++, i++)
  {
    _c[2*i] = IloBoolArray(_env, m);
    _c[2*i+1] = IloBoolArray(_env, m);
    for (int j=0; j<m; j++)
    {
      _c[2*i][j] = (it->getC0() >> m-1-j) & 1;
      _c[2*i+1][j] = (it->getC1() >> m-1-j) & 1;
    }
    _homozygous[i]=it->isHomozygous();
  }
}

void IlpSolver::initIdeotype(bool swapIdeotype)
{
  // initialize ideotype, store bits in _t1 and t2, assumed _t1=_t2
  const Genotype& ideotype = _pData->getIdeotype();
  const int m = _pData->getNumberOfLoci();
  const int n = _pData->getParents().size();

  // two targets - used only in multi-target case
  for (int j=0; j<m; j++)
  {
    _t1[j] = (ideotype.getC0() >> m-1-j) & 1;
    _t2[j] = (ideotype.getC1() >> m-1-j) & 1;
  }

  // fix target bits
  for (int j=0; j<m; j++)
    _model.add(_p[(2*_options._bound)-2][j]== (swapIdeotype ? _t2[j] : _t1[j]));

  if (!ideotype.isHomozygous())
  {
    for (int j=0; j<m; j++)
      _model.add(_p[(2*_options._bound)-1][j]== (swapIdeotype ? _t1[j] : _t2[j]));
  }
}

void IlpSolver::initChromosomesFromGenotypes()
{
  // initialize x variables and corresponding constraints
  const int n = _pData->getParents().size();
  const int m = _pData->getNumberOfLoci();
  const int g = _options._fixedGen;
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();

  // backbone chromosomes are the last 2*(g-1) chromosomes
  // + first two chromosomes

  // an x variable is defined as _x[k][i],
  // _x[k][i] = 1 iff chromosome k is obtained from node i
  // if k is backbone chromosome then 0 <= i < n+k/2
  // otherwise 0 <= i < n+k/2
  // OR _options._bound - _options._fixedGen <= i < _options._bound - 3

  for (int i = 0; i < _options._bound; i++)
  {
    // add offset if g is set and i is a non-backbone node
    int offset = getNrInnerPred(i);

    _x[2*i] = IloBoolVarArray(_env, n+offset);
    _xx[2*i] = IloBoolVarArray(_env, n+offset);
    if (i < _options._bound - 1 || !homozygousIdeotype)
    {
      _x[2*i+1] = IloBoolVarArray(_env, n+offset);
      _xx[2*i+1] = IloBoolVarArray(_env, n+offset);
    }

    for (int j = 0; j < n + offset; j++)
    {
      std::stringstream ss;
      ss << "x_" << 2*i << "_" << j;
      _x[2*i][j] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_x[2*i][j]);

      if (i < _options._bound - 1 || !homozygousIdeotype)
      {
        std::stringstream sss;
        sss << "xx_" << 2*i << "_" << j;
        _xx[2*i][j] = IloBoolVar(_env, sss.str().c_str());
        _allVar.add(_xx[2*i][j]);

        std::stringstream ss2;
        ss2 << "x_" << 2*i+1 << "_" << j;
        _x[2*i+1][j] = IloBoolVar(_env, ss2.str().c_str());
        _allVar.add(_x[2*i+1][j]);

        std::stringstream sss2;
        sss2 << "xx_" << 2*i+1 << "_" << j;
        _xx[2*i+1][j] = IloBoolVar(_env, sss2.str().c_str());
        _allVar.add(_xx[2*i+1][j]);
      }
    }
  }

  IloExpr sumofxx(_env);
  for (int i = 0; i < _options._bound; i++)
  {
    int offset = getNrInnerPred(i);
    for (int j = 0; j < n + offset; j++)
    {
      if (i < _options._bound - 1 || !homozygousIdeotype)
      {
        _model.add(_xx[2*i][j] <= _x[2*i+1][j]);
        _model.add(_xx[2*i+1][j] <= _x[2*i][j]);

        sumofxx += _xx[2*i][j];
        sumofxx -= _xx[2*i+1][j];
      }
    }
    _model.add(sumofxx == 0);
    sumofxx.clear();
  }

  // every chromosome has exactly one predecessor
  IloExpr sumofinedges1(_env);
  IloExpr sumofinedges2(_env);
  for (int i=0; i < _options._bound; i++)
  {
    int offset = getNrInnerPred(i);
    for (int j=0; j < n + offset; j++)
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
    for (int i = 0; i < _options._bound; i++)
    {
      // backbone nodes don't matter here, their out-degree is by definition >= 1
      if (!isBackboneNode(i))
      {
        // TODO: maybe upper bound for j can be tighter... DONE!
        // (REMEMBER: last node is only connected to the last but one)
        const int chromosomeUB = homozygousIdeotype ? 2 * _options._bound - 1 : 2 * _options._bound;
        for (int j = 2*(i+1); j < chromosomeUB; j++)
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
    for (int h = 0; h < _options._fixedGen - 2; h++)
    {
      int i = _options._bound - h - 1;
      // we skip n+i-1, because we connect to n+i-1
      for (int j = 0; j < n+i-1; j++)
      {
        _model.add(_x[2*i][j] == 0);
      }
      _model.add(_x[2*i][n+i-1] == 1);
    }

    // connection from first backbone node to 2nd backbone node
    int i = _options._bound - _options._fixedGen + 1;
    for (int j = 0; j < n + i; j++)
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
  const int n = _pData->getParents().size();
  const int m = _pData->getNumberOfLoci();
  const int g = _options._fixedGen;
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();

  // y-Vars, encodes bit origins
  // y[i][j] = { 0, if j-th bit of chromosome i originates 
  //                from the upper chromosome of the node defined by x
  //             1, if j-th bit of chromosome i originates
  //                from the lower chromosome of the node defined by x
  const int chromosomeUB = homozygousIdeotype ? 2 * _options._bound - 1 : 2 * _options._bound;
  for (int i = 0; i < chromosomeUB; i++)
  {
    _y[i] = IloBoolVarArray(_env, m);
    for (int j = 0; j < m; j++)
    {
      std::stringstream ss;
      ss << "y_" << i << "_" << j;
      _y[i][j] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_y[i][j]);
    }
  }

  // the bit variable
  for (int i = 0; i < chromosomeUB; i++)
  {
    _p[i] = IloBoolVarArray(_env, m);
    for (int t = 0; t < m; t++)
    {
      std::stringstream ss;
      ss << "p_" << i << "_" << t;
      _p[i][t] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_p[i][t]);
    }
  }

  // is locus heterozygous, 
  // (not defined for the target node, as it is homozygous by definition)
  const int genotypeUB = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  for (int i = 0; i < genotypeUB; i++)
  {
    _bt[i] = IloBoolVarArray(_env, m);
    for (int t = 0; t < m; t++)
    {
      std::stringstream ss;
      ss << "bt_" << i << "_" << t;
      _bt[i][t] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_bt[i][t]);
    }
  }

  for (int i = 0; i < genotypeUB; i++)
  {
    for (int j = 0; j < m; j++)
    {
      // TODO reicht das? keine <= Ungleichung notwendig?
      // Ja das reicht! Being homozygous is beneficial
      _model.add(_bt[i][j] >= _p[2*i][j] - _p[2*i+1][j]);
      _model.add(_bt[i][j] >= _p[2*i+1][j] - _p[2*i][j]);
    }
  }

  // is node heterozygous
  for (int i = 0; i < genotypeUB; i++)
  {
    std::stringstream ss;
    ss << "het_" << i;
    _het[i] = IloBoolVar(_env, ss.str().c_str());
    _allVar.add(_het[i]);
    for (int j = 0; j < m; j++)
      _model.add(_het[i] >= _bt[i][j]);
  }

  // selfing lemma:
  // if a node i is homozygous (i.e. _het[i]==0), 
  // then it's two predecessors must be the same
  for (int i = 0; i < genotypeUB; i++)
  {
    int offset = getNrInnerPred(i);
    for (int j = 0; j < n + offset; j++)
    {
      _model.add(_x[2*i][j] - _x[2*i+1][j] <= _het[i]);
      _model.add(_x[2*i+1][j] - _x[2*i][j] <= _het[i]);
    }
  }

  // _input_bit[i][j][t] = 1 iff 1-allele at locus j from chr i comes from chr t
  for (int i = 0; i < chromosomeUB; i++)
  {
    _input_bit[i] = BoolVarMatrix(_env, m);

    int offset = 2*getNrInnerPred(i/2);

    for (int j = 0; j < m; j++)
    {
      _input_bit[i][j] = IloBoolVarArray(_env, 2*n + offset);
      for (int t = 0; t < 2*n + offset; t++)
      {
        std::stringstream ss;
        ss << "input_bit_";
        ss << i << "_" << j << "_" << t;
        _input_bit[i][j][t] = IloBoolVar(_env, ss.str().c_str());
        _allVar.add(_input_bit[i][j][t]);
      }
    }
  }

  IloExpr sumof1dep(_env);
  for (int i = 0; i < chromosomeUB; i++)
  {
    for (int j = 0; j < m; j++)
    {
      // inner nodes
      int offset = 2*getNrInnerPred(i/2);
      for (int l = 0; l < offset; l++)
      {
        _model.add(_input_bit[i][j][2*n+l] <= _x[i][n + l/2]);
        _model.add(_input_bit[i][j][2*n+l] <= _p[getAbsChromosome(i/2, l)][j]);
        if ((l%2) == 0)
        {
          _model.add(_input_bit[i][j][2*n+l] <= 1-_y[i][j]);
          _model.add(_input_bit[i][j][2*n+l] >= _x[i][l/2+n] + _p[getAbsChromosome(i/2,l)][j] + (1-_y[i][j])-2);
        }
        else
        {
          _model.add(_input_bit[i][j][2*n+l] <= _y[i][j]);
          _model.add(_input_bit[i][j][2*n+l] >= _x[i][l/2+n] + _p[getAbsChromosome(i/2,l)][j] + _y[i][j]-2);
        }
      }

      // leaves
      for (int l = 0; l < 2*n; l++)
      {
        _model.add(_input_bit[i][j][l] <= _c[l][j] * _x[i][l/2]);
        if ((l%2) == 0)
        {
          _model.add(_input_bit[i][j][l] <= _c[l][j] * (1-_y[i][j]));
          _model.add(_input_bit[i][j][l] >= _c[l][j] * (_x[i][l/2]+(1-_y[i][j])-1));
        }
        else
        {
          _model.add(_input_bit[i][j][l] <= _c[l][j] * _y[i][j]);
          _model.add(_input_bit[i][j][l] >= _c[l][j] * (_x[i][l/2]+_y[i][j]-1));
        }
      }  
      for (int l=0; l < 2*n + offset; l++)
      {
        sumof1dep += _input_bit[i][j][l];
      }
      _model.add(sumof1dep == _p[i][j]);
      sumof1dep.clear();
    }
  }
  sumof1dep.end();
}

void IlpSolver::initUsefulCross()
{ 
  const int m = _pData->getNumberOfLoci();
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const int chromosomeUB = homozygousIdeotype ? 2 * _options._bound - 1 : 2 * _options._bound;

  IloExpr sumofcollected(_env);
  for (int i=0; i < chromosomeUB; i++)
  {
    for (int j=0; j<m-2; j++)
    {
      for (int l=j+1; l<m-1; l++)
      {
        for (int t=j+1; t<l+1; t++)
          sumofcollected+=_p[i][t];
        _model.add(_d[i][j]+_d[i][l]-sumofcollected<=1);
        sumofcollected.clear();
      }
    }

    for (int j=0; j<m-1; j++)
    {
      for (int t=0; t<j+1; t++)
        sumofcollected+=_p[i][t];
      _model.add(_d[i][j]-sumofcollected<=0);
      sumofcollected.clear();
    }      
    for (int j=0; j<m-1; j++)
    {
      for (int t=j+1; t<m; t++)
        sumofcollected+=_p[i][t];
      _model.add(_d[i][j]-sumofcollected<=0);
      sumofcollected.clear();
    }
  }
  sumofcollected.end();
}

void IlpSolver::initDiffParents()
{
  // TODO: fix this!
  const int n = _pData->getParents().size();
  const int g = _options._fixedGen;
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const int genotypeUB = homozygousIdeotype ? _options._bound - 1 : _options._bound;

  for (int i=0; i < genotypeUB; i++)
  {
    for (int j=0; j<n; j++)
      _model.add(_x[2*i][j]+_x[2*i+1][j]<=1);

    int offset = g <= 1 || i >= _options._bound-g ? 0 : g-2;
    for (int j=0; j<i+offset; j++)
      _model.add(_x[2*i][j+n]+_x[2*i+1][j+n]<=1);
  }
}

void IlpSolver::initTree()
{
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
      if (_t1[p] == _c[2*i][p] || _t1[p] == _c[2*i+1][p])
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
        if (cover[p] == 1 && _c[2*j][p] == _t1[p])
        {
          sum1+=_x[2*i][j];
        }
        if (i < _options._bound-1 && cover[p] == 1 && _c[2*j+1][p] == _t1[p])
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
        _model.add(_p[2*i][p] == sum1);
        if (i < _options._bound-1)
        {
          _model.add(_p[2*i+1][p] == sum2);
          // homozygous => 0,0
          _model.add(_p[2*i][p] + _p[2*i+1][p] <= _bt[i][p]);
        }
      }
      else
      {
        // if locus p at chromosome 2i is unique, we must keep it
        _model.add(_p[2*i][p] >= sum1);
        if (i < _options._bound-1)
          // if locus p at chromosome 2i+1 is unique, we must keep it
          _model.add(_p[2*i+1][p] >= sum2);
      }
      sum1.clear();
      sum2.clear();
    }
  }
  sum1.end();
  sum2.end();
}

void IlpSolver::initObligatoryLeaves()
{
  const int m = _pData->getNumberOfLoci();
  const int n = _pData->getParents().size();
  const GenotypeSet& parents = _pData->getParents();
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();

  // outdegree of obligatory leaves is at least one
  bool* obligatory = new bool[n];
  for (int i=0; i<n; i++)
    obligatory[i]=false;
 
  for (int j=0; j<m; j++)
  {
    int i=0;
    bool unique_covered = false;
    int obl_parent = -1;
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


  const int chromosomeUB = homozygousIdeotype ? 2 * _options._bound - 1 : 2 * _options._bound;
  IloExpr sumofoutedges(_env);
  for (int i=0; i<n; i++)
  {
    if (!obligatory[i])
      continue;
    for (int j=0; j < chromosomeUB; j++)
    {    
      sumofoutedges+=_x[j][i];
    }
    _model.add(sumofoutedges>=1);       
    sumofoutedges.clear();
  }    
  sumofoutedges.end();
  delete obligatory;
}

void IlpSolver::initObjectiveGen()
{
  const int n = _pData->getParents().size();
  const int bound = _options._bound;
  const int g = _options._fixedGen;

  _obj += _pData->getCostGen() * _r[bound-1];
  for (int i = 0; i < bound; i++)
  {
    std::stringstream ss;
    ss << "r_" << i;
    _r[i].setName(ss.str().c_str());
    _allVar.add(_r[i]);

    int offset = getNrInnerPred(i);
    for (int j = 0; j < offset; j++)
    {
      _model.add(_r[i] >= _r[getAbsNode(i,j)] + 1 - (g*(1-_x[2*i][n+j])));
      if (i < bound-1)
        _model.add(_r[i] >= _r[getAbsNode(i,j)] + 1 - (g*(1-_x[2*i+1][n+j])));
    }
  }
}

void IlpSolver::initObjectiveCross()
{
  _obj += _pData->getCostNode() * _options._bound;
}

void IlpSolver::initObjectivePop()
{
  const int n = _pData->getParents().size();
  const int m = _pData->getNumberOfLoci();
  const DoubleMatrix& RM = _pData->getRM();
  const int g = _options._fixedGen;
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();

  // product var for heterozygous and x
  // _hetx[k][i] = _het[i]*_x[k][i]
  for (int i = 0; i < _options._bound; i++)
  {
    int offset = getNrInnerPred(i);
    _hetx[2*i] = IloBoolVarArray(_env, n + offset);
    if (i <_options._bound - 1 || !homozygousIdeotype)
      _hetx[2*i+1] = IloBoolVarArray(_env, n + offset);

    // leaves
    for (int j = 0; j < n; j++)
    { 
      std::stringstream ss;
      ss << "hetx_" << 2*i << "_" << j;
      _hetx[2*i][j] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_hetx[2*i][j]);
      if (i < _options._bound - 1 || !homozygousIdeotype)
      {
        std::stringstream ss;
        ss << "hetx_" << 2*i + 1 << "_" << j;
        _hetx[2*i+1][j] = IloBoolVar(_env, ss.str().c_str());
        _allVar.add(_hetx[2*i+1][j]);
      }
      if (_homozygous[j])
      {
        _model.add(_hetx[2*i][j]==0);
        if (i<_options._bound-1 || !homozygousIdeotype)
          _model.add(_hetx[2*i+1][j]==0);
      }
      else
      {
        _hetx[2*i][j] = _x[2*i][j];
        if (i<_options._bound-1 || !homozygousIdeotype)
          _hetx[2*i+1][j] = _x[2*i+1][j];
      }
    }

    // inner nodes
    for (int j = n; j < n + offset; j++)
    {
      std::stringstream ss;
      ss << "hetx_" << 2*i << "_" << j;
      _hetx[2*i][j] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_hetx[2*i][j]);
      if (i < _options._bound - 1 || !homozygousIdeotype)
      {
        std::stringstream ss;
        ss << "hetx_" << 2*i + 1 << "_" << j;
        _hetx[2*i+1][j] = IloBoolVar(_env, ss.str().c_str());
        _allVar.add(_hetx[2*i+1][j]);
      }

      _model.add(_hetx[2*i][j]<=_het[getAbsNode(i,j-n)]);
      _model.add(_hetx[2*i][j]<=_x[2*i][j]);
      _model.add(_hetx[2*i][j]>=_het[getAbsNode(i,j-n)]+_x[2*i][j]-1);

      if (i < _options._bound - 1 || !homozygousIdeotype)
      {
        _model.add(_hetx[2*i+1][j]<=_het[getAbsNode(i,j-n)]);
        _model.add(_hetx[2*i+1][j]<=_x[2*i+1][j]);
        _model.add(_hetx[2*i+1][j]>=_het[getAbsNode(i,j-n)]+_x[2*i+1][j]-1);
      }
    }
  }

  // did jth bit in string i result from a crossover event?
  const int chromosomeUB = homozygousIdeotype ? 2 * _options._bound - 1 : 2 * _options._bound;
  for (int i = 0; i < chromosomeUB; i++)
  {
    _d[i] = IloBoolVarArray(_env, m-1);
    for (int j = 0; j < m - 1; j++)
    {
      std::stringstream ss;
      ss << "d_" << i << "_" << j;
      _d[i][j] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_d[i][j]);
    }
  }
  for (int i = 0; i < chromosomeUB; i++)
  {
    for (int j = 1; j < m; j++)
    {
      _model.add(_d[i][j-1] >= _y[i][j]-_y[i][j-1]);
      _model.add(_d[i][j-1] >= _y[i][j-1]-_y[i][j]);
    }
  }

  // h_{i,p,q} = 1 iff loci p and q of node i are heterozygous
  // and everything in between is homozygous
  const int genotypeUB = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  IloExpr sumofr(_env);
  for (int i = 0; i < genotypeUB; i++)
  {
    _h[i] = BoolVarMatrix(_env, m - 1);
    for (int p = 0; p < m - 1; p++)
    {
      _h[i][p] = IloBoolVarArray(_env, m-p-1);
      for (int q = p + 1; q < m; q++)
      {
        std::stringstream ss;
        ss << "h_" << i << "_" << p << "_" << q-p-1;
        _h[i][p][q-p-1] = IloBoolVar(_env, ss.str().c_str());
        _allVar.add(_h[i][p][q-p-1]);
        _model.add(_h[i][p][q-p-1]<=_bt[i][p]);
        _model.add(_h[i][p][q-p-1]<=_bt[i][q]);
        for (int r = p + 1; r < q; r++)
        {
          _model.add(_h[i][p][q-p-1]<=1-_bt[i][r]);
          sumofr+=1-_bt[i][r];
        }
        _model.add(_h[i][p][q-p-1]>=_bt[i][p]+_bt[i][q]+sumofr-(q-p));
        sumofr.clear();
      }
    }
  }
  sumofr.end();

  // product var of _h and _x:
  // _hx[k][i][p][q] = _h[i][p][q]*_x[k][i]
  for (int k = 0; k < chromosomeUB; k++)
  {
    int offset = getNrInnerPred(k/2);
    _hx[k] = BoolVar3Matrix(_env, n + offset);
    
    // leaves
    for (int i = 0; i < n; i++)
    {
      _hx[k][i] = BoolVarMatrix(_env, m-1);
      for (int p = 0; p < m-1; p++)
      {
        _hx[k][i][p] = IloBoolVarArray(_env, m-p-1);
        for (int q = p+1; q < m; q++)
        {
          std::stringstream ss;
          ss << "hx_" << k << "_" << i << "_" << p << "_" << q-p-1;
          _hx[k][i][p][q-p-1] = IloBoolVar(_env, ss.str().c_str());
          _allVar.add(_hx[k][i][p][q-p-1]);
          if (_c[2*i][p] == _c[2*i+1][p] ||
              _c[2*i][q] == _c[2*i+1][q] ||
              !isHomozygousBlock(i,p,q))
          {
            // _h[i][p][q] == 0
            _model.add(_hx[k][i][p][q-p-1] == 0);
          }
          else
          {
            _model.add(_hx[k][i][p][q-p-1] == _x[k][i]);
          }
        }
      }
    }

    // inner nodes
    for (int i = n; i < n + offset; i++)
    {
      _hx[k][i] = BoolVarMatrix(_env, m-1);
      for (int p = 0; p < m-1; p++)
      {
        _hx[k][i][p] = IloBoolVarArray(_env, m-p-1);
        for (int q = p+1; q < m; q++)
        {
          std::stringstream ss;
          ss << "hx_" << k << "_" << i << "_" << p << "_" << q-p-1;
          _hx[k][i][p][q-p-1] = IloBoolVar(_env, ss.str().c_str());
          _allVar.add(_hx[k][i][p][q-p-1]);
          _model.add(_hx[k][i][p][q-p-1] <= _h[getAbsNode(k/2,i-n)][p][q-p-1]);
          _model.add(_hx[k][i][p][q-p-1] <= _x[k][i]);
          _model.add(_hx[k][i][p][q-p-1] >= _h[getAbsNode(k/2,i-n)][p][q-p-1] + _x[k][i] - 1);
        }
      }
    }
  }

  // var e
  IloExpr sumofr2(_env);
  for (int k=0; k < chromosomeUB; k++)
  {
    _e[k] = IntVarMatrix(_env, m-1);
    for (int p = 0; p < m-1; p++)
    {
      _e[k][p] = IloIntVarArray(_env, m-p-1);
      for (int q = p+1; q < m; q++)
      {
        std::stringstream ss;
        ss << "e_" << k << "_" << p << "_" << q-p-1;
        _e[k][p][q-p-1] = IloIntVar(_env, 0, m-1, ss.str().c_str());
        _allVar.add(_e[k][p][q-p-1]);
        for (int r = p + 1; r <= q; r++)
        {
          sumofr2 += _d[k][r-1]; // mind the r-1
        }
        _model.add(_e[k][p][q-p-1] == sumofr2);
        sumofr2.clear();
      }
    }
  }
  sumofr2.end();

  // product var of _h, _x and _e:
  // _hxe[k][i][p][q] = _h[i][p][q]*_x[k][i]*_e[k][p][q]
  for (int k=0; k < chromosomeUB; k++)
  {
    int offset = getNrInnerPred(k/2);
    _hxe[k] = BoolVar3Matrix(_env, n + offset);
    
    for (int i = 0; i < n; i++)
    {
      _hxe[k][i] = BoolVarMatrix(_env, m - 1);
      for (int p = 0; p < m - 1; p++)
      {
        _hxe[k][i][p] = IloBoolVarArray(_env, m-p-1);
        for (int q = p+1; q < m; q++)
        {
          std::stringstream ss;
          ss << "hxe_" << k << "_" << i << "_" << p << "_" << q-p-1;
          _hxe[k][i][p][q-p-1] = IloBoolVar(_env, ss.str().c_str());
          _allVar.add(_hxe[k][i][p][q-p-1]);
          if (_c[2*i][p] == _c[2*i+1][p] ||
              _c[2*i][q] == _c[2*i+1][q] ||
              !isHomozygousBlock(i,p,q))
          {
            // _h[i][p][q] == 0
            _model.add(_hxe[k][i][p][q-p-1] == 0);
          }
          else
          {
            _model.add(_hxe[k][i][p][q-p-1] <= _x[k][i]);
            _model.add(_hxe[k][i][p][q-p-1] <= _e[k][p][q-p-1]);
            _model.add(_hxe[k][i][p][q-p-1] >= _x[k][i] + _e[k][p][q-p-1] - 1);
          }
        }
      }
    }
    for (int i = n; i < n + offset; i++)
    {
      _hxe[k][i] = BoolVarMatrix(_env, m-1);
      for (int p = 0; p < m-1; p++)
      {
        _hxe[k][i][p] = IloBoolVarArray(_env, m-p-1);
        for (int q = p+1; q < m; q++)
        {
          std::stringstream ss;
          ss << "hxe_" << k << "_" << i << "_" << p << "_" << q-p-1;
          _hxe[k][i][p][q-p-1] = IloBoolVar(_env, ss.str().c_str());
          _allVar.add(_hxe[k][i][p][q-p-1]);
          // TODO, maybe better in terms of h and x
          _model.add(_hxe[k][i][p][q-p-1] <= _hx[k][i][p][q-p-1]);
          _model.add(_hxe[k][i][p][q-p-1] <= _e[k][p][q-p-1]);
          _model.add(_hxe[k][i][p][q-p-1] >= _hx[k][i][p][q-p-1] + _e[k][p][q-p-1] - 1);
        }
      }
    }
  }

  // separable z var
  IloExpr sumoflogprobs(_env);
  for (int j = 0; j < _options._bound; j++)
  {
    std::stringstream ss;
    ss << "zb_" << j;
    _zb[j].setName(ss.str().c_str());
    _allVar.add(_zb[j]);

    int offset = getNrInnerPred(j);
    for (int i = 0; i < n + offset; i++)
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
        sumoflogprobs += log(0.5)*_hetx[2*j][i];
        if (j < _options._bound - 1 || !homozygousIdeotype)
          sumoflogprobs+=log(0.5)*_hetx[2*j+1][i];
      }

      for (int p = 0; p < m - 1; p++)
      {
        for (int q = p + 1; q < m; q++)
        {
          double r_pq = RM[m-p-1][m-q-1];
          sumoflogprobs += log(1-r_pq) * _hx[2*j][i][p][q-p-1];
          sumoflogprobs += log(r_pq/(1.-r_pq)) * _hxe[2*j][i][p][q-p-1];

          if (j < _options._bound - 1 || !homozygousIdeotype)
          {
            sumoflogprobs += log(1-r_pq) * _hx[2*j+1][i][p][q-p-1];
            sumoflogprobs += log(r_pq/(1.-r_pq)) * _hxe[2*j+1][i][p][q-p-1];
          }
        }
      }
    }
    if (j < _options._bound - 1 || !homozygousIdeotype)
    {
      _model.add(_zb[j] == sumoflogprobs);
      if (_options._boundNmax)
        _model.add(sumoflogprobs >= log(_pData->getProbLowerBound()));
    }
    else if (homozygousIdeotype)
    {
      // last crossing is a selfing (assumption: homozygous ideotype)
      _model.add(_zb[j] == 2 * sumoflogprobs);
      if (_options._boundNmax)
        _model.add(2*sumoflogprobs >= log(_pData->getProbLowerBound()));
    }
    sumoflogprobs.clear();
  }
  sumoflogprobs.end();

  IloExpr sumoflambdas(_env);
  for (int i=0; i < _options._bound; i++)
  {
    _lambda[i] = IloNumVarArray(_env, _breakpoint.size(), 0., 1., IloNumVar::Float);
    _lambda[i][_breakpoint.size()-1] = IloBoolVar(_env);
    std::stringstream ss;
    ss << "lambda_" << i << "_" << _breakpoint.size() - 1;
    _lambda[i][_breakpoint.size()-1].setName(ss.str().c_str());
    _allVar.add(_lambda[i][_breakpoint.size()-1]);

    for (int j = 0; j < _breakpoint.size()-1; j++)
    {
      std::stringstream ss;
      ss << "lambda_" << i << "_" << j;
      _lambda[i][j].setName(ss.str().c_str());
      _allVar.add(_lambda[i][j]);
      sumoflambdas+=_lambda[i][j];
    }
    _model.add(sumoflambdas == 1 - _lambda[i][_breakpoint.size()-1]);
    sumoflambdas.clear();

    for (int j = 0; j < _breakpoint.size(); j++)
    {
      sumoflambdas+=_breakpoint[j]*_lambda[i][j];
    }
    _model.add(_zb[i]==sumoflambdas);
    sumoflambdas.clear();
  }
  sumoflambdas.end();

  IloExpr sumofpopsize(_env);
  for (int i = 0; i < _options._bound; i++)
  {
    for (int j = 0; j < _breakpoint.size(); j++)
    {
      if (j == _breakpoint.size() - 1)
      {
        assert(_breakpoint[j] == 0 && _pData->getGamma() >= .5);
        sumofpopsize+=_lambda[i][j];
      }
      else
      {
        sumofpopsize+=(log(1-_pData->getGamma())/log(1-exp(_breakpoint[j])))*_lambda[i][j];
      }
    }
  }

  _obj += _pData->getCostCrossover() * sumofpopsize;
  sumofpopsize.clear();
  sumofpopsize.end();
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
  for (int i = 0; i < _options._bound-1; i++)
  {
    genotypeVector.push_back(new Genotype(parseGenotype(i)));
  }
  genotypeVector.push_back(new Genotype(_pData->getIdeotype()));

  // reconstruct arcs and nodes
  bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  for (int i = 0; i < _options._bound; i++)
  {
    if (!genotypeVector[n+i])
      continue;

    int p1 = -1, p2 = -1;
    int offset = getNrInnerPred(i);
    for (int j = 0; j < n + offset; j++)
    {
      if (_pCplex->getValue(_x[2*i][j]))
      {
        assert(p1 == -1);
        if (j>=n+i)
          p1 = n+getAbsNode(i, j-n);
        else
          p1 = j;
      }
      if (homozygousIdeotype)
      {
        if (i != _options._bound - 1 && _pCplex->getValue(_x[2*i+1][j]))
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
        if (_pCplex->getValue(_x[2*i+1][j]))
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

void IlpSolver::printEdges(Node node, 
                           BoolNodeMap& visited, 
                           std::ostream& out) const
{
  static const int nLoci = _pData->getNumberOfLoci();

  if (!visited[node])
  {
    for (InArcIt a(_G, node); a != lemon::INVALID; ++a)
    {
      Node parent = _G.source(a);
      out << '\t';
      _genotype[parent].printGenotype(nLoci, false, out, "");
      out << " -> ";
      _genotype[node].printGenotype(nLoci, false, out, "");
      out << " [label=\"" 
        << (_options._printProb ? exp(parseLogProb(_idx[node])) : _pop[node]) 
        << "\"];" << std::endl;
      printEdges(parent, visited, out);
    }
    visited[node] = true;
  }
}
