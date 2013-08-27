/*
 * ilpsolver.h
 *
 *  Created on: 26-apr-2011
 *      Author: M. El-Kebir
 */

#ifndef ILPSOLVER_H_
#define ILPSOLVER_H_

#include <ilcplex/ilocplex.h>
#include <limits>
#include "common/data.h"
#include "crossingschedule.h"

class IlpSolver : public CrossingSchedule
{
public:
  struct Options {
    bool _tree;
    bool _noSelfing;
    bool _usefulCross;
    bool _verbose;
    size_t _bound;     // fixes the number of internal nodes
    size_t _fixedGen;  // fixes the number of generations
    bool _boundNmax;
    int _nof_segs;
    bool _servin;
    double _upperBoundObj;
    double _upperBoundPop;

    /// Default constructor with default values
    Options()
      : _tree(false)
      , _noSelfing(false)
      , _usefulCross(false)
      , _verbose(false)
      , _bound(1)
      , _fixedGen(1)
      , _boundNmax(false)
      , _nof_segs(7)
      , _servin(false)
      , _upperBoundObj(std::numeric_limits<double>::max())
      , _upperBoundPop(std::numeric_limits<double>::max())
    {
    }
  };

  typedef enum {
    CSO_SOLVER_OPTIMAL,
    CSO_SOLVER_INFEASIBLE,
    CSO_SOLVER_TIME_LIMIT_INFEASIBLE,
    CSO_SOLVER_TIME_LIMIT_FEASIBLE
  } SolverStatus;

protected:
  typedef IloArray<IloBoolVarArray> BoolVarMatrix;
  typedef IloArray<BoolVarMatrix> BoolVar3Matrix;
  typedef IloArray<BoolVar3Matrix> BoolVar4Matrix;
  typedef IloArray<IloNumVarArray> NumVarMatrix;
  typedef IloArray<NumVarMatrix> NumVar3Matrix;
  typedef IloArray<IloIntVarArray> IntVarMatrix;
  typedef IloArray<IntVarMatrix> IntVar3Matrix;

  const Options& _options;
  IdxMap _idx;
  IloEnv _env;
  IloModel _model;
  IloNumVarArray _allVar;
  IloExpr _obj;
  IloCplex* _pCplex;
  /// Piece-wise linear function breakpoints: log probablities
  DoubleVector _breakpoint;
  /// Piece-wise linear function breakpoints: population sizes
  DoubleVector _N;
  /// C^*_1
  IloBoolArray _cs1;
  /// C^*_2
  IloBoolArray _cs2;
  /// Encodes parental genotypes
  IloArray<IloBoolArray> _c;
  /// Encodes whether a parental genotype is homozygous
  IloBoolArray _homozygous;
  /// Encodes whether a chromosome originates from a certain node
  BoolVarMatrix _x;
  /// Encodes g_{k,p,l} where l corresponds to a parental chromosome
  BoolVar3Matrix _g;
  /// y[i][j] = { 0, if j-th bit of chromosome i originates
  ///                from the upper chromosome of the node defined by x
  ///             1, if j-th bit of chromosome i originates
  ///                from the lower chromosome of the node defined by x
  BoolVarMatrix _y;
  /// Allele value at given chromosome and locus
  BoolVarMatrix _a;
  /// Encodes whether a locus at given node is heterozygous
  BoolVarMatrix _at;
  /// Encodes whether a given node is heterozygous
  IloBoolVarArray _h;
  /// Encodes product variable for h and x
  BoolVarMatrix _hx;
  /// Encodes whether allele is the result of a crossover event
  BoolVarMatrix _d;
  /// Probability for obtaining a node
  IloNumVarArray _p;
  /// Lambda for piece-wise linear approximation of pop size
  NumVarMatrix _lambda;
  /// Encodes depth of a node
  IloIntVarArray _r;
  /// Encodes b_{i,p,q}
  BoolVar3Matrix _b;
  /// Encodes product of b and x
  BoolVar4Matrix _bx;
  /// Encodes product of b and x and z
  BoolVar4Matrix _bxz;
  /// Encodes z_{k,p,q}:
  /// the alleles at loci p and q of chromosome k are the result of a crossover
  IntVar3Matrix _z;
  /// Uniquely covered
  BoolVarMatrix _uc;
  /// Product variable uc*x
  BoolVar3Matrix _ucx;

  virtual void initSegments();
  virtual void initParentalGenotypes();
  virtual void initIdeotype(bool swapIdeotype);
  virtual void initChromosomesFromGenotypes();
  virtual void initAllelesFromChromosomes();
  virtual void initUsefulCross();
  virtual void initNoSelfing();
  virtual void initTree();
  virtual void initObligatoryLeaves();
  virtual void initGen();
  virtual void initPop();
  virtual void initBounds();
  virtual void initObj();
  virtual void initPopExpr(IloExpr& expr) const;

  bool isHomozygousBlock(size_t i, size_t p, size_t q) const;
  void printA() const;
  void printAt() const;
  void printB() const;
  void printZ() const;
  void printX() const;
  void printY() const;
  void printBX() const;
  void printBXZ() const;
  void printG() const;
  void printP() const;
  virtual void printLambda() const;
  void printInnerNodes() const;

  Genotype parseGenotype(size_t i) const;
  virtual double parsePop(size_t i) const;
  virtual double parseProb1(size_t i) const;
  virtual double parseProb2(size_t i) const;
  unsigned long parseGen(size_t i) const;
  void constructDAG();
  size_t getRelNode(size_t base_i, size_t backbone_i) const;
  size_t getAbsChromosome(size_t base_i, size_t backbone_k) const;
  size_t getAbsNode(size_t base_i, size_t backbone_i) const;

  bool isBackboneNode(size_t i) const;
  size_t getNrInnerPred(size_t i) const;

public:
  IlpSolver(const Data* pData, const Options& options);
  virtual ~IlpSolver();
  void init(bool swapIdeotype = false);
  virtual SolverStatus solve(bool feasibility, int timeLimit);
  double getObjectiveValue() const;
};

inline bool IlpSolver::isBackboneNode(size_t i) const
{
  return i == 0 ||
      (_options._bound - _options._fixedGen < i && i < _options._bound);
}

inline size_t IlpSolver::getNrInnerPred(size_t i) const
{
  if (isBackboneNode(i))
  {
    return i;
  }
  else
  {
    // -3, because the last and second-last backbone nodes
    // as well as first backbone node cannot be predecessors
    return i + _options._fixedGen - 3;
  }
}

inline size_t IlpSolver::getAbsChromosome(size_t i, size_t k) const
{
  if (k < 2*i)
    return k;
  else
    return k - 2*i + 2*(_options._bound - _options._fixedGen + 1);
}

inline size_t IlpSolver::getAbsNode(size_t i, size_t j) const
{
  if (j < i)
    return j;
  else
    return j - i + _options._bound - _options._fixedGen + 1;
}

inline size_t IlpSolver::getRelNode(size_t base_i, size_t backbone_i) const
{
  abort();
  // TODO
  const size_t firstBackboneNode = _options._bound - _options._fixedGen;
  return base_i + backbone_i - firstBackboneNode;
}

inline double IlpSolver::getObjectiveValue() const
{
  return _pCplex->getObjValue();
}

inline void IlpSolver::printInnerNodes() const
{
  const size_t m = _pData->getNumberOfLoci();
  for (size_t i = 0; i < _options._bound; i++)
  {
    std::cout << "// _r[" << i << "] = "
      << _pCplex->getValue(_r[i]) << ": ";

    if (i < _options._bound - 1)
      parseGenotype(i).printGenotype(m);
    else
      _pData->getIdeotype().printGenotype(m);
  }
}

inline void IlpSolver::printA() const
{
  const size_t chromosomeUB = _pData->isIdeotypeHomozygous() ? 2 * _options._bound - 1 : 2 * _options._bound;
  const size_t m = _pData->getNumberOfLoci();
  for (size_t k = 0; k < chromosomeUB; k++)
  {
    for (size_t p = 0; p < m; p++)
    {
      std::cout << "// a[" << k << "][" << p << "] = "
                << _pCplex->getValue(_a[k][p]) << std::endl;
    }
  }
}

inline void IlpSolver::printAt() const
{
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const size_t genotypeUB = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  const size_t m = _pData->getNumberOfLoci();
  for (size_t i = 0; i < genotypeUB; i++)
  {
    for (size_t p = 0; p < m; p++)
    {
      std::cout << "// at[" << i << "][" << p << "] = "
                << _pCplex->getValue(_at[i][p]) << std::endl;
    }
  }
}

inline void IlpSolver::printB() const
{
  const bool homozygousIdeotype = _pData->isIdeotypeHomozygous();
  const size_t genotypeUB = homozygousIdeotype ? _options._bound - 1 : _options._bound;
  const size_t m = _pData->getNumberOfLoci();
  for (size_t i = 0; i < genotypeUB; i++)
  {
    for (size_t p = 0; p < m-1; p++)
    {
      for (size_t q = p+1; q < m; q++)
      {
        std::cout << "// b[" << i << "][" << p << "][" << q << "] = "
          << _pCplex->getValue(_b[i][p][q-p-1]) << std::endl;
      }
    }
  }
}

inline void IlpSolver::printP() const
{
  for (size_t i = 0; i < _options._bound; i++)
  {
    double val = _pCplex->getValue(_p[i]);
    std::cout << "// p[" << i << "] = "
      << val << " " << exp(val) << std::endl;
  }
}

inline void IlpSolver::printLambda() const
{
  for (size_t i = 0; i < _options._bound; i++)
  {
    for (size_t j = 0; j < _breakpoint.size(); j++)
    {
      double val = _pCplex->getValue(_lambda[i][j]);
      std::cout << "// lambda[" << i << "][" << j
                << "] = " << val << std::endl;
    }
  }
}

inline void IlpSolver::printZ() const
{
  const size_t m = _pData->getNumberOfLoci();
  const size_t chromosomeUB = _pData->isIdeotypeHomozygous() ? 2 * _options._bound - 1 : 2 * _options._bound;
  for (size_t k = 0; k < chromosomeUB; k++)
  {
    for (size_t p = 0; p < m-1; p++)
    {
      for (size_t q = p+1; q < m; q++)
      {
        std::cout << "// z[" << k << "][" << p << "][" << q << "] = "
          << _pCplex->getValue(_z[k][p][q-p-1]) << std::endl;
      }
    }
  }
}

inline void IlpSolver::printY() const
{
  const size_t chromosomeUB = _pData->isIdeotypeHomozygous() ? 2 * _options._bound - 1 : 2 * _options._bound;
  const size_t m = _pData->getNumberOfLoci();
  for (size_t k = 0; k < chromosomeUB; k++)
  {
    for (size_t p = 0; p < m; p++)
    {
      std::cout << "// y[" << k << "][" << p << "] = "
        << _pCplex->getValue(_y[k][p]) << std::endl;
    }
  }
}

inline void IlpSolver::printX() const
{
  const size_t chromosomeUB = _pData->isIdeotypeHomozygous() ? 2 * _options._bound - 1 : 2 * _options._bound;
  const size_t n = _pData->getParents().size();
  for (size_t k = 0; k < chromosomeUB; k++)
  {
    for (size_t i = 0; i < n + getNrInnerPred(k/2); i++)
    {
      std::cout << "// x[" << k << "][" << i << "] = "
        << _pCplex->getValue(_x[k][i]) << std::endl;
    }
  }
}

inline void IlpSolver::printG() const
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
        std::cout << "// g[" << k << "][" << p << "][" << i << "] = "
          << _pCplex->getValue(_g[k][p][i]) << std::endl;
      }
    }
  }
}

inline void IlpSolver::printBX() const
{
  const size_t chromosomeUB = _pData->isIdeotypeHomozygous() ? 2 * _options._bound - 1 : 2 * _options._bound;
  const size_t m = _pData->getNumberOfLoci();
  const size_t  n = _pData->getParents().size();
  for (size_t k = 0; k < chromosomeUB; k++)
  {
    for (size_t i = 0; i < n+k/2; i++)
    {
      for (size_t p = 0; p < m-1; p++)
      {
        for (size_t q = p+1; q < m; q++)
        {
          int val1 = _pCplex->getValue(_bx[k][i][p][q-p-1]);
          int val2 = 0;

          if (i < n)
          {
            if (_c[2*i][p] == _c[2*i+1][p] || 
                _c[2*i][q] == _c[2*i+1][q] ||
                !isHomozygousBlock(i,p,q))
            {
              // _h[i][p][q] == 0
              val2 = 0;
            }
            else
            {
              val2 = _pCplex->getValue(_x[k][i]);
            }
          }
          else
          {
            val2 = _pCplex->getValue(_b[i-n][p][q-p-1]) * _pCplex->getValue(_x[k][i]);
          }

          std::cout << "// bx[" << k << "][" << i << "]["
            << p << "][" << q << "] = " << val1 
            << (val1 == val2 ? "=" : "!=") << val2 << std::endl;
        }
      }
    }
  }
}

inline void IlpSolver::printBXZ() const
{
  const size_t chromosomeUB = _pData->isIdeotypeHomozygous() ? 2 * _options._bound - 1 : 2 * _options._bound;
  const size_t m = _pData->getNumberOfLoci();
  const size_t  n = _pData->getParents().size();
  for (size_t k = 0; k < chromosomeUB; k++)
  {
    for (size_t i = 0; i < n+k/2; i++)
    {
      for (size_t p = 0; p < m-1; p++)
      {
        for (size_t q = p+1; q < m; q++)
        {
          int val1 = _pCplex->getValue(_bxz[k][i][p][q-p-1]);
          int val2 = 0;

          if (i < n)
          {
            if (_c[2*i][p] == _c[2*i+1][p] || 
                _c[2*i][q] == _c[2*i+1][q] ||
                !isHomozygousBlock(i,p,q))
            {
              // _h[i][p][q] == 0
              val1 = 0;
            }
            else
            {
              val2 = _pCplex->getValue(_x[k][i]) * _pCplex->getValue(_z[k][p][q-p-1]);
            }
          }
          else
          {
            val2 = _pCplex->getValue(_b[i-n][p][q-p-1]) *
              _pCplex->getValue(_x[k][i]) * _pCplex->getValue(_z[k][p][q-p-1]);
          }
          std::cout << "// bxe[" << k << "][" << i << "]["
            << p << "][" << q << "] = " << val1 
            << (val1 == val2 ? "=" : "!=") << val2 << std::endl;
        }
      }
    }
  }
}

inline bool IlpSolver::isHomozygousBlock(size_t i, size_t p, size_t q) const
{
  const size_t m = _pData->getNumberOfLoci();
  const size_t n = _pData->getParents().size();

  assert(0 <= i && i < n);
  assert(0 <= p && p < q && q < m);

  for (size_t r = p + 1; r < q; r++)
  {
    if (_c[2*i][r] != _c[2*i+1][r])
      return false;
  }
  return true;
}

inline Genotype IlpSolver::parseGenotype(size_t i) const
{
  assert(0 <= i && i < _options._bound);
  const size_t m = _pData->getNumberOfLoci();

  int c0 = 0, c1 = 0;
  for (size_t l = 0; l < m; l++)
  {
    bool b0 = _pCplex->getValue(_a[2*i][l]) != 0;
    bool b1 = _pCplex->getValue(_a[2*i+1][l]) != 0;

    c0 |= (int)b0 << (m - l - 1);
    c1 |= (int)b1 << (m - l - 1);
  }

  return Genotype(c0, c1);
}

inline double IlpSolver::parsePop(size_t i) const
{
  assert(0 <= i && i < _options._bound);

  double pop = 0;
  for (size_t j = 0; j < _breakpoint.size(); j++)
  {
    pop += _pCplex->getValue(_lambda[i][j]) * _N[j];
  }

  return pop;
}

inline double IlpSolver::parseProb1(size_t i) const
{
  assert(0 <= i && i < _options._bound);
  return exp(_pCplex->getValue(_p[i]));
}

inline double IlpSolver::parseProb2(size_t i) const
{
  assert(0 <= i && i < _options._bound);
  return 0;
}

inline unsigned long IlpSolver::parseGen(size_t i) const
{
  assert(0 <= i && i < _options._bound);
  return lround(_pCplex->getValue(_r[i]));
}

#endif /* ILPSOLVER_H_ */
