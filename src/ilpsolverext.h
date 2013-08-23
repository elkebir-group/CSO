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
  /// Piece-wise linear function breakpoints (pop): log probablities
  DoubleVector _breakpoint2;
  /// Mu for piece-wise linear approximation of probablities (_p)
  NumVarMatrix _mu;
  /// Nu for piece-wise linear approximation of probablities (_pp)
  NumVarMatrix _nu;
  /// Decision variables
  BoolVarMatrix _b_mu;
  BoolVarMatrix _b_nu;

  virtual void initPopExpr(IloExpr& expr) const;
  virtual void initSegments();
  virtual void initChromosomesFromGenotypes();
  virtual void initAllelesFromChromosomes();
  virtual void initPop();
  virtual void initObj();
  void printH() const;
  void printF() const;
  void printXX() const;
  void printPP() const;
  void printGG() const;
  void printMu() const;
  void printNu() const;
  virtual double parseProb2(int i) const;

public:
  IlpSolverExt(const Data* pData, const Options& options);
  virtual ~IlpSolverExt();
  virtual SolverStatus solve(bool feasibility, int timeLimit);
};

inline void IlpSolverExt::printXX() const
{
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const int genotypeUB = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  const int n = _pData->getParents().size();
  for (int j = 0; j < genotypeUB; j++)
  {
    for (int k = 2*j; k <= 2*j+1; k++)
    {
      for (int i = 0; i < n + getNrInnerPred(j); i++)
      {
        std::cout << "// xx[" << k << "][" << i << "] = "
          << _pCplex->getValue(_xx[k][i]) << std::endl;
      }
    }
  }
}

inline void IlpSolverExt::printMu() const
{
  for (size_t i = 0; i < _options._bound; i++)
  {
    for (size_t j = 0; j < _breakpoint.size(); j++)
    {
      double val = _pCplex->getValue(_mu[i][j]);
      std::cout << "// mu[" << i << "][" << j
                << "] = " << val << std::endl;
    }
  }
}

inline void IlpSolverExt::printNu() const
{
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const size_t genotypeUB2 = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  for (size_t i = 0; i < genotypeUB2; i++)
  {
    for (size_t j = 0; j < _breakpoint.size(); j++)
    {
      double val = _pCplex->getValue(_nu[i][j]);
      std::cout << "// nu[" << i << "][" << j
                << "] = " << val << std::endl;
    }
  }
}

inline void IlpSolverExt::printH() const
{
  for (int i = 0; i < _options._bound; i++)
  {
    std::cout << "// h[" << i << "] = "
              << _pCplex->getValue(_h[i]) << std::endl;
  }
}

inline void IlpSolverExt::printF() const
{
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const int genotypeUB2 = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  for (int i = 0; i < genotypeUB2; i++)
  {
    std::cout << "// f[" << i << "] = "
              << _pCplex->getValue(_f[i]) << std::endl;
  }
}


inline void IlpSolverExt::printPP() const
{
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const int genotypeUB2 = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  for (int i = 0; i < genotypeUB2; i++)
  {
    double val = _pCplex->getValue(_pp[i]);
    std::cout << "// pp[" << i << "] = "
      << val << " " << exp(val) << std::endl;
  }
}

inline void IlpSolverExt::printGG() const
{
  const int chromosomeUB = _pData->isIdeotypeHomozygous() ? 2 * _options._bound - 1 : 2 * _options._bound;
  const int n = _pData->getParents().size();
  const int m = _pData->getNumberOfLoci();
  for (int k = 0; k < chromosomeUB; k++)
  {
    for (int p = 0; p < m; p++)
    {
      for (int i = 0; i < 2*(n + getNrInnerPred(k/2)); i++)
      {
        std::cout << "// gg[" << k << "][" << p << "][" << i << "] = "
          << _pCplex->getValue(_gg[k][p][i]) << std::endl;
      }
    }
  }
}

inline double IlpSolverExt::parseProb2(int i) const
{
  assert(0 <= i && i < _options._bound);

  if (_pData->isIdeotypeHomozygous() && i == _options._bound - 1)
 {
    return 0;
  }

  //if (_pCplex->getValue(_f[i]))
  //{
  //  return exp(_pCplex->getValue(_pp[i]));
  //}
  //else
  {
    return 0;
  }
}


#endif // ILPSOLVEREXT_H
