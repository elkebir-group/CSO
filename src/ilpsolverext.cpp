/*
 * ilpsolverext.cpp
 *
 *  Created on: 22-aug-2013
 *      Author: M. El-Kebir
 */

#include "ilpsolverext.h"

IlpSolverExt::IlpSolverExt(const Data* pData, const Options& options)
  : IlpSolver(pData, options)
  , _xx(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-2 : 2 * options._bound)
  , _f(_env, pData->isIdeotypeHomozygous() ? options._bound - 1 : options._bound)
  , _gg(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-2 : 2 * options._bound)
  , _yy(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-2 : 2 * options._bound)
  , _hxx(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-2 : 2 * options._bound)
  , _dd(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-2 : 2 * options._bound)
  , _pp(_env, pData->isIdeotypeHomozygous() ? options._bound - 1 : options._bound, -IloInfinity, 0.0, IloNumVar::Float)
  , _bxx(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-2 : 2 * options._bound)
  , _bxxzz(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-2 : 2 * options._bound)
  , _zz(_env, pData->isIdeotypeHomozygous() ? (2 * options._bound)-2 : 2 * options._bound)
  , _mu(_env, options._bound)
  , _nu(_env, options._bound)
  , _b_lambda1(_env, options._bound)
  , _b_lambda2(_env, options._bound)
  , _breakpoint1()
  , _breakpoint2()
  , _N()
  , _lambda(_env, options._bound)
{
}

IlpSolverExt::~IlpSolverExt()
{
  _xx.end();
  _gg.end();
  _yy.end();
  _hxx.end();
  _dd.end();
  _pp.end();
  _bxx.end();
  _bxxzz.end();
  _zz.end();
  _mu.end();
  _nu.end();
  _b_lambda1.end();
  _b_lambda2.end();
  _lambda.end();
}

void IlpSolverExt::initSegments()
{
  _breakpoint1.clear();
  _breakpoint2.clear();
  double boundary = log(1-pow(1-_pData->getGamma(),1./_pData->getPopMax()));
  double seglength = (log(.5)-boundary)/_options._nof_segs;

  _breakpoint2.push_back(0);
  for (int i = 0; i <= _options._nof_segs; i++)
  {
    double val = boundary + i * seglength;
    _breakpoint1.push_back(val);
    _breakpoint2.push_back(val);
  }
  _breakpoint1.push_back(0);

  _N.clear();
  for (size_t j = 0; j < _breakpoint1.size(); j++)
  {
    const double p1 = exp(_breakpoint1[j]);

    const size_t k_max = j == _breakpoint1.size() - 1 ? 0 : j + 1;
    _N.push_back(DoubleVector(k_max + 1, 0));
    for (size_t k = 0; k <= k_max; k++)
    {
      const double p2 = k == 0 ? 0 : exp(_breakpoint2[k]);
      const double pop = std::max(1., log(1-_pData->getGamma())/log(1-(p1+p2)));
      _N[j][k] = pop;
    }
  }
}

void IlpSolverExt::initPopExpr(IloExpr& expr) const
{
  expr.clear();
  for (size_t i = 0; i < _options._bound; i++)
  {
    for (size_t j = 0; j < _breakpoint1.size(); j++)
    {
      const size_t k_max = j == _breakpoint1.size() - 1 ? 0 : j + 1;
      for (size_t k = 0; k <= k_max; k++)
      {
        expr += _lambda[i][j][k] * _N[j][k];
      }
    }
  }
}

IlpSolver::SolverStatus IlpSolverExt::solve(bool feasibility, int timeLimit)
{
  SolverStatus stat = IlpSolver::solve(feasibility, timeLimit);
  if (stat != CSO_SOLVER_INFEASIBLE &&
      stat != CSO_SOLVER_TIME_LIMIT_INFEASIBLE &&
      _options._verbose)
  {
    //std::cout << "//" << std::endl;
    //printH();
    //std::cout << "//" << std::endl;
    //printF();
    //std::cout << "//" << std::endl;
    //printXX();
    //std::cout << "//" << std::endl;
    //printX();
    //std::cout << "//" << std::endl;
    //printYY();
    //std::cout << "//" << std::endl;
    //printGG();
    //std::cout << "//" << std::endl;
    //printPP();
    //std::cout << "//" << std::endl;
    //printBLambda2();
    //std::cout << "//" << std::endl;
  }
  return stat;
}

void IlpSolverExt::initObj()
{
  IloExpr sumoflambdas(_env);
  IloExpr sumoflambdas_only_pp(_env);
  IloExpr sumoflambdas_p(_env);
  IloExpr sumoflambdas_pp(_env);

  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const size_t genotypeUB2 = homozygousIdeotype ? _options._bound - 1 : _options._bound;

  for (size_t i = 0; i < _options._bound; i++)
  {
    _lambda[i] = NumVarMatrix(_env, _breakpoint1.size());

    for (size_t j = 0; j < _breakpoint1.size(); j++)
    {
      const size_t k_max = j == _breakpoint1.size() - 1 ? 0 : j + 1;

      _lambda[i][j] = j == _breakpoint1.size() - 1 ?
            IloNumVarArray(_env, k_max + 1, 0, 1, IloNumVar::Int) :
            IloNumVarArray(_env, k_max + 1, 0., 1., IloNumVar::Float);
      for (size_t k = 0; k <= k_max; k++)
      {
        std::stringstream ss;
        ss << "lambda_" << i << "_" << j << "_" << k;
        _lambda[i][j][k].setName(ss.str().c_str());
        _allVar.add(_lambda[i][j][k]);

        sumoflambdas_p += _breakpoint1[j] * _lambda[i][j][k];
        sumoflambdas_pp += _breakpoint2[k] * _lambda[i][j][k];
        sumoflambdas += _lambda[i][j][k];

        if (k != 0)
        {
          sumoflambdas_only_pp += _lambda[i][j][k];
        }
      }
    }

    _model.add(sumoflambdas == 1);

    if (i < genotypeUB2 - 1 || !homozygousIdeotype)
    {
      _model.add(sumoflambdas_only_pp == _f[i]);
      _model.add(sumoflambdas_pp == _pp[i]);
    }
    else if (homozygousIdeotype)
    {
      _model.add(sumoflambdas_only_pp == 0);
    }

    _model.add(sumoflambdas_p == _p[i]);

    sumoflambdas.clear();
    sumoflambdas_only_pp.clear();
    sumoflambdas_p.clear();
    sumoflambdas_pp.clear();
  }

  for (size_t i = 0; i < _options._bound; i++)
  {
    //_b_lambda1[i] = IloBoolVarArray(_env, _breakpoint1.size() - 1);
    _b_lambda2[i] = IloBoolVarArray(_env, _breakpoint2.size() - 1);

    //IloExpr sum_b_lambda1(_env);
    IloExpr sum_b_lambda2(_env);
    //for (size_t j = 0; j < _breakpoint1.size() - 1; j++)
    //{
    //  std::stringstream ss;
    //  ss << "b_lambda1_" << i << "_" << j;
    //  _b_lambda1[i][j] = IloBoolVar(_env);
    //  _b_lambda1[i][j].setName(ss.str().c_str());
    //  _allVar.add(_b_lambda1[i][j]);
//
    //  sum_b_lambda1 += _b_lambda1[i][j];
    //}

    for (size_t k = 0; k < _breakpoint2.size() - 1; k++)
    {
      std::stringstream ss2;
      ss2 << "b_lambda2_" << i << "_" << k;
      _b_lambda2[i][k] = IloBoolVar(_env);
      _b_lambda2[i][k].setName(ss2.str().c_str());
      _allVar.add(_b_lambda2[i][k]);

      sum_b_lambda2 += _b_lambda2[i][k];
    }

    //for (size_t j = 0; j < _breakpoint1.size(); j++)
    //{
    //  const size_t k_max = j == _breakpoint1.size() - 1 ? 0 : j + 1;
    //  for (size_t k = 0; k <= k_max; k++)
    //  {
    //    sumoflambdas += _lambda[i][j][k];
    //  }
//
    //  if (j == 0)
    //  {
    //    _model.add(sumoflambdas <= _b_lambda1[i][0]);
    //  }
    //  else if (j != _breakpoint1.size() - 1)
    //  {
    //    _model.add(sumoflambdas <= _b_lambda1[i][j-1] + _b_lambda1[i][j]);
    //  }
    //  else
    //  {
    //    _model.add(sumoflambdas <= _b_lambda1[i][j-1]);
    //  }
    //  sumoflambdas.clear();
    //}
    //_model.add(sum_b_lambda1 == 1);

    for (size_t k = 0; k < _breakpoint2.size(); k++)
    {
      const size_t j_min = k == 0 ? 0 : k - 1;
      const size_t j_max = k == 0 ? _breakpoint1.size() - 1 : _breakpoint1.size() - 2;
      for (size_t j = j_min; j <= j_max; j++)
      {
        sumoflambdas += _lambda[i][j][k];
      }

      if (k == 0)
      {
        _model.add(sumoflambdas <= _b_lambda2[i][0]);
      }
      else if (k != _breakpoint2.size() - 1)
      {
        _model.add(sumoflambdas <= _b_lambda2[i][k-1] + _b_lambda2[i][k]);
      }
      else
      {
        _model.add(sumoflambdas <= _b_lambda2[i][k-1]);
      }
      sumoflambdas.clear();
    }
    _model.add(sum_b_lambda2 == 1);
  }

  initPopExpr(_obj);
}

void IlpSolverExt::initChromosomesFromGenotypes()
{
  IlpSolver::initChromosomesFromGenotypes();

  // initialize xx variables and corresponding constraints
  const size_t n = _pData->getParents().size();
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const size_t genotypeUB2 = homozygousIdeotype ? _options._bound - 1 : _options._bound;

  for (size_t i = 0; i < genotypeUB2; i++)
  {
    // add offset if g is set and i is a non-backbone node
    int offset = getNrInnerPred(i);

    _xx[2*i] = IloBoolVarArray(_env, n+offset);
    _xx[2*i+1] = IloBoolVarArray(_env, n+offset);

    for (size_t j = 0; j < n + offset; j++)
    {
      std::stringstream sss;
      sss << "xx_" << 2*i << "_" << j;
      _xx[2*i][j] = IloBoolVar(_env, sss.str().c_str());
      _allVar.add(_xx[2*i][j]);

      std::stringstream sss2;
      sss2 << "xx_" << 2*i+1 << "_" << j;
      _xx[2*i+1][j] = IloBoolVar(_env, sss2.str().c_str());
      _allVar.add(_xx[2*i+1][j]);
    }
  }

  // \sum_{j=1}^{i-1} xx[2i][j] == \sum_{j=1}^{i-1} xx[2i+1][j]
  // _xx[2*i][j] <= _x[2*i+1][j]
  // _xx[2*i+1][j] <= _x[2*i][j]
  IloExpr sumofxx(_env);
  for (size_t i = 0; i < _options._bound; i++)
  {
    size_t offset = getNrInnerPred(i);
    for (size_t j = 0; j < n + offset; j++)
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

  // \sum_{j=1}^{i-1} xx[2i][j] == f[i]
  for (size_t i = 0; i < genotypeUB2; i++)
  {
    std::stringstream ss;
    ss << "f_" << i;
    _f[i] = IloBoolVar(_env, ss.str().c_str());
    _allVar.add(_f[i]);

    size_t offset = getNrInnerPred(i);
    for (size_t j = 0; j < n + offset; j++)\
    {
      sumofxx += _xx[2*i][j];
    }
    _model.add(sumofxx == _f[i]);
    sumofxx.clear();
  }
}

void IlpSolverExt::initAllelesFromChromosomes()
{
  IlpSolver::initAllelesFromChromosomes();

  const size_t n = _pData->getParents().size();
  const size_t m = _pData->getNumberOfLoci();
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();

  const size_t genotypeUB2 = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  const size_t chromosomeUB2 = homozygousIdeotype ? 2 * _options._bound - 2 : 2 * _options._bound;

  // _xx[2*i][j] <= _h[i]
  // _xx[2*i+1][j] <= _h[i]
  for (size_t i = 0; i < genotypeUB2; i++)
  {
    size_t offset = getNrInnerPred(i);
    for (size_t j = 0; j < n + offset; j++)
    {
      _model.add(_xx[2*i][j] <= _h[i]);
      _model.add(_xx[2*i+1][j] <= _h[i]);
    }
  }

  // y-Vars, encodes bit origins
  // yy[i][j] = { 0, if j-th bit of chromosome i originates
  //                 from the upper chromosome of the node defined by xx
  //              1, if j-th bit of chromosome i originates
  //                 from the lower chromosome of the node defined by xx
  for (size_t i = 0; i < chromosomeUB2; i++)
  {
    _yy[i] = IloBoolVarArray(_env, m);

    for (size_t j = 0; j < m; j++)
    {
      std::stringstream sss;
      sss << "yy_" << i << "_" << j;
      _yy[i][j] = IloBoolVar(_env, sss.str().c_str());
      _allVar.add(_yy[i][j]);
    }
  }

  // _gg[i][j][t] = 1 iff 1-allele at locus j from chr i comes from chr t
  for (size_t i = 0; i < chromosomeUB2; i++)
  {
    _gg[i] = BoolVarMatrix(_env, m);

    size_t offset = 2*getNrInnerPred(i/2);
    for (size_t j = 0; j < m; j++)
    {
      _gg[i][j] = IloBoolVarArray(_env, 2*n + offset);

      for (size_t t = 0; t < 2*n + offset; t++)
      {
        std::stringstream sss;
        sss << "gg_";
        sss << i << "_" << j << "_" << t;
        _gg[i][j][t] = IloBoolVar(_env, sss.str().c_str());
        _allVar.add(_gg[i][j][t]);
      }
    }
  }

  for (size_t i = 0; i < chromosomeUB2; i++)
  {
    for (size_t j = 0; j < m; j++)
    {
      // inner nodes
      size_t offset = 2*getNrInnerPred(i/2);
      for (size_t l = 0; l < offset; l++)
      {
        _model.add(_gg[i][j][2*n+l] <= _xx[i][n + l/2]);
        _model.add(_gg[i][j][2*n+l] <= _a[getAbsChromosome(i/2, l)][j]);
        if ((l%2) == 0)
        {
          _model.add(_gg[i][j][2*n+l] <= 1-_yy[i][j]);
          _model.add(_gg[i][j][2*n+l] >= _xx[i][l/2+n] + _a[getAbsChromosome(i/2,l)][j] + (1-_yy[i][j])-2);

          _model.add(_gg[i][j][2*n+l] + _gg[i][j][2*n+l+1] <= _a[i][j]);
          _model.add(_gg[i][j][2*n+l] + _gg[i][j][2*n+l+1] <= _xx[i][l/2+n]);
          _model.add(_gg[i][j][2*n+l] + _gg[i][j][2*n+l+1] >= _a[i][j] + _xx[i][l/2+n] - 1);
        }
        else
        {
          _model.add(_gg[i][j][2*n+l] <= _yy[i][j]);
          _model.add(_gg[i][j][2*n+l] >= _xx[i][l/2+n] + _a[getAbsChromosome(i/2,l)][j] + _yy[i][j]-2);
        }
      }

      // leaves
      for (size_t l = 0; l < 2*n; l++)
      {
        _model.add(_gg[i][j][l] <= _c[l][j] * _xx[i][l/2]);
        if ((l%2) == 0)
        {
          _model.add(_gg[i][j][l] <= _c[l][j] * (1-_yy[i][j]));
          _model.add(_gg[i][j][l] >= _c[l][j] * (_xx[i][l/2]+(1-_yy[i][j])-1));

          _model.add(_gg[i][j][l] + _gg[i][j][l+1] <= _a[i][j]);
          _model.add(_gg[i][j][l] + _gg[i][j][l+1] <= _xx[i][l/2]);
          _model.add(_gg[i][j][l] + _gg[i][j][l+1] >= _a[i][j] + _xx[i][l/2] - 1);
        }
        else
        {
          _model.add(_gg[i][j][l] <= _c[l][j] * _yy[i][j]);
          _model.add(_gg[i][j][l] >= _c[l][j] * (_xx[i][l/2]+_yy[i][j]-1));
        }
      }
    }
  }
}

void IlpSolverExt::initPop()
{
  IlpSolver::initPop();

  const size_t n = _pData->getParents().size();
  const size_t m = _pData->getNumberOfLoci();
  const DoubleMatrix& RM = _pData->getRM();
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();

  const size_t genotypeUB2 = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  const size_t chromosomeUB2 = homozygousIdeotype ? _options._bound * 2 - 2 : _options._bound * 2;

  // product var for heterozygous and x
  // _hxx[k][i] = _hxx[i]*_xx[k][i]
  for (size_t i = 0; i < genotypeUB2; i++)
  {
    size_t offset = getNrInnerPred(i);
    _hxx[2*i] = IloBoolVarArray(_env, n + offset);
    _hxx[2*i+1] = IloBoolVarArray(_env, n + offset);

    // leaves
    for (size_t j = 0; j < n; j++)
    {
      std::stringstream sss;
      sss << "hxx_" << 2*i << "_" << j;
      _hxx[2*i][j] = IloBoolVar(_env, sss.str().c_str());
      _allVar.add(_hxx[2*i][j]);

      std::stringstream sss2;
      sss2 << "hxx_" << 2*i + 1 << "_" << j;
      _hxx[2*i+1][j] = IloBoolVar(_env, sss2.str().c_str());
      _allVar.add(_hxx[2*i+1][j]);

      if (_homozygous[j])
      {
        _model.add(_hxx[2*i][j]==0);
        _model.add(_hxx[2*i+1][j]==0);
      }
      else
      {
        _hxx[2*i][j] = _xx[2*i][j];
        _hxx[2*i+1][j] = _xx[2*i+1][j];
      }
    }

    // inner nodes
    for (size_t j = n; j < n + offset; j++)
    {
      std::stringstream sss;
      sss << "hxx_" << 2*i << "_" << j;
      _hxx[2*i][j] = IloBoolVar(_env, sss.str().c_str());
      _allVar.add(_hxx[2*i][j]);

      std::stringstream sss2;
      sss2 << "hxx_" << 2*i + 1 << "_" << j;
      _hxx[2*i+1][j] = IloBoolVar(_env, sss2.str().c_str());
      _allVar.add(_hxx[2*i+1][j]);

      _model.add(_hxx[2*i][j]<=_h[getAbsNode(i,j-n)]);
      _model.add(_hxx[2*i][j]<=_xx[2*i][j]);
      _model.add(_hxx[2*i][j]>=_h[getAbsNode(i,j-n)]+_xx[2*i][j]-1);

      _model.add(_hxx[2*i+1][j]<=_h[getAbsNode(i,j-n)]);
      _model.add(_hxx[2*i+1][j]<=_xx[2*i+1][j]);
      _model.add(_hxx[2*i+1][j]>=_h[getAbsNode(i,j-n)]+_xx[2*i+1][j]-1);
    }
  }

  // did jth bit in string i result from a crossover event?
  for (size_t i = 0; i < chromosomeUB2; i++)
  {
    _dd[i] = IloBoolVarArray(_env, m-1);
    for (size_t j = 0; j < m - 1; j++)
    {
      std::stringstream ss;
      ss << "dd_" << i << "_" << j;
      _dd[i][j] = IloBoolVar(_env, ss.str().c_str());
      _allVar.add(_dd[i][j]);
    }
  }

  for (size_t i = 0; i < chromosomeUB2; i++)
  {
    for (size_t j = 1; j < m; j++)
    {
      _model.add(_dd[i][j-1] >= _yy[i][j]-_yy[i][j-1]);
      _model.add(_dd[i][j-1] >= _yy[i][j-1]-_yy[i][j]);
    }
  }

  // product var of _b and _xx:
  // _bxx[k][i][p][q] = _b[i][p][q]*_xx[k][i]
  for (size_t k = 0; k < chromosomeUB2; k++)
  {
    size_t offset = getNrInnerPred(k/2);
    _bxx[k] = BoolVar3Matrix(_env, n + offset);

    // leaves
    for (size_t i = 0; i < n; i++)
    {
      _bxx[k][i] = BoolVarMatrix(_env, m-1);
      for (size_t p = 0; p < m-1; p++)
      {
        _bxx[k][i][p] = IloBoolVarArray(_env, m-p-1);
        for (size_t q = p+1; q < m; q++)
        {
          std::stringstream ss;
          ss << "bxx_" << k << "_" << i << "_" << p << "_" << q-p-1;
          _bxx[k][i][p][q-p-1] = IloBoolVar(_env, ss.str().c_str());
          _allVar.add(_bxx[k][i][p][q-p-1]);
          if (_c[2*i][p] == _c[2*i+1][p] ||
              _c[2*i][q] == _c[2*i+1][q] ||
              !isHomozygousBlock(i,p,q))
          {
            // _bxx[i][p][q] == 0
            _model.add(_bxx[k][i][p][q-p-1] == 0);
          }
          else
          {
            _model.add(_bxx[k][i][p][q-p-1] == _xx[k][i]);
          }
        }
      }
    }

    // inner nodes
    for (size_t i = n; i < n + offset; i++)
    {
      _bxx[k][i] = BoolVarMatrix(_env, m-1);
      for (size_t p = 0; p < m-1; p++)
      {
        _bxx[k][i][p] = IloBoolVarArray(_env, m-p-1);
        for (size_t q = p+1; q < m; q++)
        {
          std::stringstream ss;
          ss << "bxx_" << k << "_" << i << "_" << p << "_" << q-p-1;
          _bx[k][i][p][q-p-1] = IloBoolVar(_env, ss.str().c_str());
          _allVar.add(_bxx[k][i][p][q-p-1]);
          _model.add(_bxx[k][i][p][q-p-1] <= _b[getAbsNode(k/2,i-n)][p][q-p-1]);
          _model.add(_bxx[k][i][p][q-p-1] <= _xx[k][i]);
          _model.add(_bxx[k][i][p][q-p-1] >= _b[getAbsNode(k/2,i-n)][p][q-p-1] + _xx[k][i] - 1);
        }
      }
    }
  }

  // var zz: zz[k][p][q] = 1 iff the alleles at loci p and q of chr k
  // are the result of a crossover
  IloExpr sumofr2(_env);
  for (size_t k = 0; k < chromosomeUB2; k++)
  {
    _zz[k] = IntVarMatrix(_env, m-1);

    for (size_t p = 0; p < m-1; p++)
    {
      _zz[k][p] = IloIntVarArray(_env, m-p-1);
      for (size_t q = p+1; q < m; q++)
      {
        std::stringstream ss;
        ss << "zz_" << k << "_" << p << "_" << q-p-1;
        _zz[k][p][q-p-1] = IloIntVar(_env, 0, m-1, ss.str().c_str());
        _allVar.add(_zz[k][p][q-p-1]);
        for (size_t r = p + 1; r <= q; r++)
        {
          sumofr2 += _dd[k][r-1]; // mind the r-1
        }
        _model.add(_zz[k][p][q-p-1] == sumofr2);
        sumofr2.clear();
      }
    }
  }
  sumofr2.end();

  // product var of _b, _xx and _zz:
  // _bxxzz[k][i][p][q] = _b[i][p][q]*_xx[k][i]*_zz[k][p][q]
  for (size_t k=0; k < chromosomeUB2; k++)
  {
    size_t offset = getNrInnerPred(k/2);
    _bxxzz[k] = IntVar3Matrix(_env, n + offset);

    for (size_t i = 0; i < n; i++)
    {
      _bxxzz[k][i] = IntVarMatrix(_env, m - 1);
      for (size_t p = 0; p < m - 1; p++)
      {
        _bxxzz[k][i][p] = IloIntVarArray(_env, m-p-1);
        for (size_t q = p+1; q < m; q++)
        {
          std::stringstream ss;
          ss << "bxxzz_" << k << "_" << i << "_" << p << "_" << q-p-1;
          _bxxzz[k][i][p][q-p-1] = IloIntVar(_env, 0, m-1, ss.str().c_str());
          _allVar.add(_bxxzz[k][i][p][q-p-1]);
          if (_c[2*i][p] == _c[2*i+1][p] ||
              _c[2*i][q] == _c[2*i+1][q] ||
              !isHomozygousBlock(i,p,q))
          {
            // _b[i][p][q] == 0
            _model.add(_bxxzz[k][i][p][q-p-1] == 0);
          }
          else
          {
            _model.add(_bxxzz[k][i][p][q-p-1] <= _xx[k][i] * static_cast<int>(m-1));
            _model.add(_bxxzz[k][i][p][q-p-1] <= _zz[k][p][q-p-1]);
            _model.add(_bxxzz[k][i][p][q-p-1] >= _xx[k][i] + (1./(m-1)) * _zz[k][p][q-p-1] - 1);
          }
        }
      }
    }
    for (size_t i = n; i < n + offset; i++)
    {
      _bxxzz[k][i] = IntVarMatrix(_env, m-1);
      for (size_t p = 0; p < m-1; p++)
      {
        _bxxzz[k][i][p] = IloIntVarArray(_env, m-p-1);
        for (size_t q = p+1; q < m; q++)
        {
          std::stringstream ss;
          ss << "bxxzz_" << k << "_" << i << "_" << p << "_" << q-p-1;
          _bxxzz[k][i][p][q-p-1] = IloIntVar(_env, 0, m-1, ss.str().c_str());
          _allVar.add(_bxxzz[k][i][p][q-p-1]);
          _model.add(_bxxzz[k][i][p][q-p-1] <= _bxx[k][i][p][q-p-1] * static_cast<int>(m-1));
          _model.add(_bxxzz[k][i][p][q-p-1] <= _zz[k][p][q-p-1]);
          _model.add(_bxxzz[k][i][p][q-p-1] >= _bxx[k][i][p][q-p-1] + (1./(m-1)) * _zz[k][p][q-p-1] - 1);
        }
      }
    }
  }

  // separable pp var
  IloExpr sumoflogprobs(_env);
  for (size_t j = 0; j < genotypeUB2; j++)
  {
    std::stringstream ss;
    ss << "pp_" << j;
    _pp[j].setName(ss.str().c_str());
    _allVar.add(_pp[j]);

    size_t offset = getNrInnerPred(j);
    for (size_t i = 0; i < n + offset; i++)
    {
      if (i < n && !_homozygous[i])
      {
        // factor of 0.5 if parent is heterozygous
        sumoflogprobs += log(0.5) * _xx[2*j][i];
        if (j < _options._bound - 1 || !homozygousIdeotype)
          sumoflogprobs += log(0.5) * _xx[2*j+1][i];
      }
      else if (i >= n)
      {
        // factor of 0.5 of parent is heterozygous
        sumoflogprobs += log(0.5)*_hxx[2*j][i];
        if (j < _options._bound - 1 || !homozygousIdeotype)
          sumoflogprobs+=log(0.5)*_hxx[2*j+1][i];
      }

      for (size_t p = 0; p < m - 1; p++)
      {
        for (size_t q = p + 1; q < m; q++)
        {
          double r_pq = RM[p][q];
          sumoflogprobs += log(1-r_pq) * _bxx[2*j][i][p][q-p-1];
          sumoflogprobs += log(r_pq/(1.-r_pq)) * _bxxzz[2*j][i][p][q-p-1];

          if (j < _options._bound - 1 || !homozygousIdeotype)
          {
            sumoflogprobs += log(1-r_pq) * _bxx[2*j+1][i][p][q-p-1];
            sumoflogprobs += log(r_pq/(1.-r_pq)) * _bxxzz[2*j+1][i][p][q-p-1];
          }
        }
      }
    }
    _model.add(_pp[j] == sumoflogprobs);
    sumoflogprobs.clear();
  }
}
