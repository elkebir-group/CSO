/*
 * ilpsolverext.h
 *
 *  Created on: 22-aug-2013
 *      Author: M. El-Kebir
 */

#ifndef ILPSOLVEREXT_H
#define ILPSOLVEREXT_H

#include "ilpsolver.h"

class IlpSolverExt : public IlpSolver
{
protected:
  /// Shadow x
  BoolVarMatrix _xx;
  /// Denotes whether node i is ambiguous
  IloBoolVarArray _f;
  /// Shadow g
  BoolVar3Matrix _gg;
  /// Shadow y
  BoolVarMatrix _yy;
  /// Shadow hx
  BoolVarMatrix _hxx;
  /// Shadow d
  BoolVarMatrix _dd;
  /// Shadow p
  IloNumVarArray _pp;
  /// Shadow bx
  BoolVar4Matrix _bxx;
  /// Shadow bxz
  BoolVar4Matrix _bxxzz;
  /// Shadow z
  IntVar3Matrix _zz;
  /// Mu for piece-wise linear approximation of probablities (_p)
  NumVarMatrix _mu;
  /// Nu for piece-wise linear approximation of probablities (_pp)
  NumVarMatrix _nu;
  /// Decision variables
  BoolVarMatrix _b_lambda1;
  BoolVarMatrix _b_lambda2;
  /// Piece-wise linear function breakpoints: log probablities
  DoubleVector _breakpoint1;
  DoubleVector _breakpoint2;
  /// Piece-wise linear function breakpoints: population sizes
  DoubleMatrix _N;
  /// Lambda for piece-wise linear approximation of pop size
  NumVar3Matrix _lambda;

  virtual void initPopExpr(IloExpr& expr) const;
  virtual void initSegments();
  virtual void initChromosomesFromGenotypes();
  virtual void initAllelesFromChromosomes();
  virtual void initPop();
  virtual void initObj();
  void printH() const;
  void printF() const;
  void printYY() const;
  void printXX() const;
  void printPP() const;
  void printGG() const;
  void printBLambda1() const;
  void printBLambda2() const;
  virtual void printLambda() const;
  virtual double parsePop(size_t i) const;
  virtual double parseProb2(size_t i) const;

public:
  IlpSolverExt(const Data* pData, const Options& options);
  virtual ~IlpSolverExt();
  virtual SolverStatus solve(bool feasibility, int timeLimit);
};

inline void IlpSolverExt::printLambda() const
{
  for (size_t i = 0; i < _options._bound; i++)
  {
    for (size_t j = 0; j < _breakpoint1.size(); j++)
    {
      const size_t k_max = j == _breakpoint1.size() - 1 ? 0 : j + 1;
      for (size_t k = 0; k <= k_max; k++)
      {
        double val = _pCplex->getValue(_lambda[i][j][k]);
        std::cout << "// lambda[" << i << "][" << j
                  << "][" << k << "] = " << val
                  << "\t" << _breakpoint1[j]
                  << "\t" << _breakpoint2[k]
                  << "\t" << _N[j][k]
                  << std::endl;
      }
    }
  }
}

inline double IlpSolverExt::parsePop(size_t i) const
{
  assert(0 <= i && i < _options._bound);

  double pop = 0;
  for (size_t j = 0; j < _breakpoint1.size(); j++)
  {
    const size_t k_max = j == _breakpoint1.size() - 1 ? 0 : j + 1;
    for (size_t k = 0; k <= k_max; k++)
    {
      pop += _pCplex->getValue(_lambda[i][j][k]) * _N[j][k];
    }
  }

  return pop;
}

inline void IlpSolverExt::printYY() const
{
  const size_t chromosomeUB = _pData->isIdeotypeHomozygous() ? 2 * _options._bound - 1 : 2 * _options._bound;
  const size_t m = _pData->getNumberOfLoci();
  for (size_t k = 0; k < chromosomeUB; k++)
  {
    for (size_t p = 0; p < m; p++)
    {
      std::cout << "// yy[" << k << "][" << p << "] = "
                << _pCplex->getValue(_yy[k][p]) << std::endl;
    }
  }
}

inline void IlpSolverExt::printXX() const
{
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const size_t genotypeUB = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  const size_t n = _pData->getParents().size();
  for (size_t j = 0; j < genotypeUB; j++)
  {
    for (size_t k = 2*j; k <= 2*j+1; k++)
    {
      for (size_t i = 0; i < n + getNrInnerPred(j); i++)
      {
        std::cout << "// xx[" << k << "][" << i << "] = "
                  << _pCplex->getValue(_xx[k][i]) << std::endl;
      }
    }
  }
}

inline void IlpSolverExt::printBLambda1() const
{
  for (size_t i = 0; i < _options._bound; i++)
  {
    for (size_t j = 0; j < _breakpoint1.size() - 1; j++)
    {
      std::cout << "// b_lambda1[" << i << "][" << j
                << "] = " << _pCplex->getValue(_b_lambda1[i][j]) << std::endl;
    }
  }
}

inline void IlpSolverExt::printBLambda2() const
{
  for (size_t i = 0; i < _options._bound; i++)
  {
    for (size_t k = 0; k < _breakpoint2.size() - 1; k++)
    {
      std::cout << "// b_lambda2[" << i << "][" << k
                << "] = " << _pCplex->getValue(_b_lambda2[i][k]) << std::endl;
    }
  }
}

inline void IlpSolverExt::printH() const
{
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const size_t genotypeUB2 = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  for (size_t i = 0; i < genotypeUB2; i++)
  {
    std::cout << "// h[" << i << "] = "
              << _pCplex->getValue(_h[i]) << std::endl;
  }
}

inline void IlpSolverExt::printF() const
{
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const size_t genotypeUB2 = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  for (size_t i = 0; i < genotypeUB2; i++)
  {
    std::cout << "// f[" << i << "] = "
              << _pCplex->getValue(_f[i]) << std::endl;
  }
}


inline void IlpSolverExt::printPP() const
{
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const size_t genotypeUB2 = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  for (size_t i = 0; i < genotypeUB2; i++)
  {
    double val = _pCplex->getValue(_pp[i]);
    std::cout << "// pp[" << i << "] = "
      << val << " " << exp(val) << std::endl;
  }
}

inline void IlpSolverExt::printGG() const
{
  const size_t chromosomeUB = _pData->isIdeotypeHomozygous() ? 2 * _options._bound - 1 : 2 * _options._bound;
  const size_t n = _pData->getParents().size();
  const size_t m = _pData->getNumberOfLoci();
  for (size_t k = 0; k < chromosomeUB; k++)
  {
    for (size_t p = 0; p < m; p++)
    {
      for (size_t i = 0; i < 2*(n + getNrInnerPred(k/2)); i++)
      {
        std::cout << "// gg[" << k << "][" << p << "][" << i << "] = "
          << _pCplex->getValue(_gg[k][p][i]) << std::endl;
      }
    }
  }
}

inline double IlpSolverExt::parseProb2(size_t i) const
{
  assert(0 <= i && i < _options._bound);

  if (_pData->isIdeotypeHomozygous() && i == _options._bound - 1)
 {
    return 0;
  }

  if (_pCplex->getValue(_f[i]))
  {
    return exp(_pCplex->getValue(_pp[i]));
  }
  else
  {
    return 0;
  }
}

#endif // ILPSOLVEREXT_H
