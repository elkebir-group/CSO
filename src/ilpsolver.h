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
    int _bound;     // fixes the number of internal nodes
    int _fixedGen;  // fixes the number of generations
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

private:
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
  /// Piece-wise linear function breakpoints
  DoubleVector _breakpoint;
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
  /// Shadow x
  BoolVarMatrix _xx;
  /// Encodes g_{k,p,l} where l corresponds to a parental chromosome
  BoolVar3Matrix _g;
  /// Shadow g
  BoolVar3Matrix _gg;
  /// y[i][j] = { 0, if j-th bit of chromosome i originates
  ///                from the upper chromosome of the node defined by x
  ///             1, if j-th bit of chromosome i originates
  ///                from the lower chromosome of the node defined by x
  BoolVarMatrix _y;
  /// Shadow y
  BoolVarMatrix _yy;
  /// Allele value at given chromosome and locus
  BoolVarMatrix _a;
  /// Encodes whether a locus at given node is heterozygous
  BoolVarMatrix _at;
  /// Encodes whether a given node is heterozygous
  IloBoolVarArray _h;
  /// Encodes product variable for h and x
  BoolVarMatrix _hx;
  /// Shadow hx
  BoolVarMatrix _hxx;
  /// Encodes whether allele is the result of a crossover event
  BoolVarMatrix _d;
  /// Shadow d
  BoolVarMatrix _dd;
  /// Probability for obtaining a node
  IloNumVarArray _p;
  /// Shadow p
  IloNumVarArray _pp;
  /// Lambda for piece-wise linear approximation of pop size
  NumVarMatrix _lambda;
  /// Encodes depth of a node
  IloIntVarArray _r;
  /// Encodes b_{i,p,q}
  BoolVar3Matrix _b;
  /// Encodes product of b and x
  BoolVar4Matrix _bx;
  /// Shadow bx
  BoolVar4Matrix _bxx;
  /// Encodes product of b and x and z
  BoolVar4Matrix _bxz;
  /// Shadow bxz
  BoolVar4Matrix _bxxzz;
  /// Encodes z_{k,p,q}:
  /// the alleles at loci p and q of chromosome k are the result of a crossover
  IntVar3Matrix _z;
  /// Shadow z
  IntVar3Matrix _zz;
  /// Uniquely covered
  BoolVarMatrix _uc;
  /// Product variable uc*x
  BoolVar3Matrix _ucx;

  void initSegments();
  void initFeasibility();
  void initObjective();
  void initParentalGenotypes();
  void initIdeotype(bool swapIdeotype);
  void initChromosomesFromGenotypes();
  void initAllelesFromChromosomes();
  void initUsefulCross();
  void initNoSelfing();
  void initTree();
  void initObligatoryLeaves();
  void initObjectiveGen();
  void initObjectiveCross();
  void initObjectivePop();
  void initBounds();
  bool isHomozygousBlock(int i, int p, int q) const;
  void printB() const;
  void printZ() const;
  void printX() const;
  void printXX() const;
  void printBX() const;
  void printBXZ() const;
  void printG() const;
  void printGG() const;
  void printP() const;
  void printPP() const;
  void printInnerNodes() const;

  Genotype parseGenotype(int i) const;
  double parsePop(int i) const;
  double parseLogProb(int i) const;
  unsigned long parseGen(int i) const;
  void constructDAG();
  int getRelNode(int base_i, int backbone_i) const;
  int getAbsChromosome(int base_i, int backbone_k) const;
  int getAbsNode(int base_i, int backbone_i) const;

  bool isBackboneNode(int i) const;
  int getNrInnerPred(int i) const;

public:
  IlpSolver(const Data* pData, const Options& options);
  ~IlpSolver();
  void init(bool swapIdeotype = false);
  SolverStatus solve(bool feasibility, int timeLimit);
  double getObjectiveValue() const;
};

inline bool IlpSolver::isBackboneNode(int i) const
{
  return i == 0 ||
      (_options._bound - _options._fixedGen < i && i < _options._bound);
}

inline int IlpSolver::getNrInnerPred(int i) const
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

inline int IlpSolver::getAbsChromosome(int i, int k) const
{
  if (k < 2*i)
    return k;
  else
    return k - 2*i + 2*(_options._bound - _options._fixedGen + 1);
}

inline int IlpSolver::getAbsNode(int i, int j) const
{
  if (j < i)
    return j;
  else
    return j - i + _options._bound - _options._fixedGen + 1;
}

inline int IlpSolver::getRelNode(int base_i, int backbone_i) const
{
  abort();
  // TODO
  const int firstBackboneNode = _options._bound - _options._fixedGen;
  return base_i + backbone_i - firstBackboneNode;
}

inline double IlpSolver::getObjectiveValue() const
{
  return _pCplex->getObjValue();
}

inline void IlpSolver::printInnerNodes() const
{
  const int m = _pData->getNumberOfLoci();
  for (int i = 0; i < _options._bound; i++)
  {
    std::cout << "// _r[" << i << "] = "
      << _pCplex->getValue(_r[i]) << ": ";

    if (i < _options._bound - 1)
      parseGenotype(i).printGenotype(m);
    else
      _pData->getIdeotype().printGenotype(m);
  }
}

inline void IlpSolver::printB() const
{
  const int m = _pData->getNumberOfLoci();
  for (int i = 0; i < _options._bound - 1; i++)
  {
    for (int p = 0; p < m-1; p++)
    {
      for (int q = p+1; q < m; q++)
      {
        std::cout << "// b[" << i << "][" << p << "][" << q << "] = "
          << _pCplex->getValue(_b[i][p][q-p-1]) << std::endl;
      }
    }
  }
}

inline void IlpSolver::printP() const
{
  for (int i = 0; i < _options._bound; i++)
  {
    double val = _pCplex->getValue(_p[i]);
    std::cout << "// p[" << i << "] = "
      << val << " " << exp(val) << std::endl;
  }
}

inline void IlpSolver::printPP() const
{
  for (int i = 0; i < _options._bound; i++)
  {
    double val = _pCplex->getValue(_pp[i]);
    std::cout << "// pp[" << i << "] = "
      << val << " " << exp(val) << std::endl;
  }
}

inline void IlpSolver::printZ() const
{
  const int m = _pData->getNumberOfLoci();
  for (int k = 0; k < 2*_options._bound - 1; k++)
  {
    for (int p = 0; p < m-1; p++)
    {
      for (int q = p+1; q < m; q++)
      {
        std::cout << "// z[" << k << "][" << p << "][" << q << "] = "
          << _pCplex->getValue(_z[k][p][q-p-1]) << std::endl;
      }
    }
  }
}

inline void IlpSolver::printXX() const
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

inline void IlpSolver::printX() const
{
  const int chromosomeUB = _pData->isIdeotypeHomozygous() ? 2 * _options._bound - 1 : 2 * _options._bound;
  const int n = _pData->getParents().size();
  for (int k = 0; k < chromosomeUB; k++)
  {
    for (int i = 0; i < n + getNrInnerPred(k/2); i++)
    {
      std::cout << "// x[" << k << "][" << i << "] = "
        << _pCplex->getValue(_x[k][i]) << std::endl;
    }
  }
}

inline void IlpSolver::printG() const
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
        std::cout << "// g[" << k << "][" << p << "][" << i << "] = "
          << _pCplex->getValue(_g[k][p][i]) << std::endl;
      }
    }
  }
}

inline void IlpSolver::printGG() const
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

inline void IlpSolver::printBX() const
{
  const int m = _pData->getNumberOfLoci();
  const int n = _pData->getParents().size();
  for (int k = 0; k < 2*_options._bound - 1; k++)
  {
    for (int i = 0; i < n+k/2; i++)
    {
      for (int p = 0; p < m-1; p++)
      {
        for (int q = p+1; q < m; q++)
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
  const int m = _pData->getNumberOfLoci();
  const int n = _pData->getParents().size();
  for (int k = 0; k < 2*_options._bound - 1; k++)
  {
    for (int i = 0; i < n+k/2; i++)
    {
      for (int p = 0; p < m-1; p++)
      {
        for (int q = p+1; q < m; q++)
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

inline bool IlpSolver::isHomozygousBlock(int i, int p, int q) const
{
  const int m = _pData->getNumberOfLoci();
  const int n = _pData->getParents().size();

  assert(0 <= i && i < n);
  assert(0 <= p && p < q && q < m);

  for (int r = p + 1; r < q; r++)
  {
    if (_c[2*i][r] != _c[2*i+1][r])
      return false;
  }
  return true;
}

inline Genotype IlpSolver::parseGenotype(int i) const
{
  assert(0 <= i && i < _options._bound);
  const int m = _pData->getNumberOfLoci();

  int c0 = 0, c1 = 0;
  for (int l = 0; l < m; l++)
  {
    bool b0 = _pCplex->getValue(_a[2*i][l]) != 0;
    bool b1 = _pCplex->getValue(_a[2*i+1][l]) != 0;

    c0 |= (int)b0 << (m - l - 1);
    c1 |= (int)b1 << (m - l - 1);
  }

  return Genotype(c0, c1);
}

inline double IlpSolver::parsePop(int i) const
{
  assert(0 <= i && i < _options._bound);

  double pop = 0;
  for(size_t j=0; j<_breakpoint.size(); j++)
  {
    if (j == _breakpoint.size() - 1)
    {
      assert(_breakpoint[j] == 0 && _pData->getGamma() >= .5);
      pop+=_pCplex->getValue(_lambda[i][j]);
    }
    else
    {
      pop+=(log(1-_pData->getGamma())/log(1-exp(_breakpoint[j])))*_pCplex->getValue(_lambda[i][j]);
    }
  }

  return pop;
}

inline double IlpSolver::parseLogProb(int i) const
{
  assert(0 <= i && i < _options._bound);

  return _pCplex->getValue(_p[i]);
}

inline unsigned long IlpSolver::parseGen(int i) const
{
  assert(0 <= i && i < _options._bound);
  return lround(_pCplex->getValue(_r[i]));
}

#endif /* ILPSOLVER_H_ */
