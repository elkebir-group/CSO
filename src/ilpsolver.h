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
    bool _diffParents;
    bool _usefulCross;
    bool _verbose;
    int _bound;     // fixes the number of internal nodes
    int _fixedGen;  // fixes the number of generations
    bool _printProb;
    bool _boundNmax;
    int _boundR;
    int _nof_segs;
    bool _servin;
    double _upperBoundObj;
    double _upperBoundPop;

    /// Default constructor with default values
    Options()
      : _tree(false)
      , _diffParents(false)
      , _usefulCross(false)
      , _verbose(false)
      , _bound(0)
      , _fixedGen(-1)
      , _printProb(false)
      , _boundNmax(false)
      , _boundR(-1)
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
  IloBoolArray _t1;
  /// C^*_2
  IloBoolArray _t2;
  /// Encodes parental genotypes
  IloArray<IloBoolArray> _c;
  /// Encodes whether a parental genotype is homozygous
  IloBoolArray _homozygous;
  /// Encodes whether a chromosome originates from a certain node
  BoolVarMatrix _x;
  /// Encodes g_{k,p,l} where l corresponds to a parental chromosome
  BoolVar3Matrix _input_bit; 
  /// y[i][j] = { 0, if j-th bit of chromosome i originates 
  ///                from the upper chromosome of the node defined by x
  ///             1, if j-th bit of chromosome i originates
  ///                from the lower chromosome of the node defined by x
  BoolVarMatrix _y;
  /// Allele value at given chromosome and locus
  BoolVarMatrix _p;
  /// Encodes whether a locus at given node is heterozygous
  BoolVarMatrix _bt;
  /// Encodes whether a given node is heterozygous
  IloBoolVarArray _het;
  /// Encodes product variable for heterozygous and x
  BoolVarMatrix _hetx;
  /// Encodes whether allele is the result of a crossover event
  BoolVarMatrix _d;
  /// Population size needed for obtaining a node
  IloNumVarArray _zb;
  /// Lambda for piece-wise linear approximation of pop size
  NumVarMatrix _lambda;
  /// Encodes depth of a node
  IloIntVarArray _r;
  /// Encodes h_{i,p,q}
  BoolVar3Matrix _h;
  /// Encodes product of h and x
  BoolVar4Matrix _hx;
  /// Encodes product of h and x and e
  BoolVar4Matrix _hxe;
  /// Encodes e_{k,p,q}: 
  /// the alleles at loci p and q of chromosome k are the result of a crossover
  IntVar3Matrix _e;
  /// Uniquely covered
  BoolVarMatrix _uc;
  /// Product variable uc*x
  BoolVar3Matrix _ucx;

  void initSegments();
  void initFeasibility();
  void initObjective();
  void initParentalGenotypes();
  void initIdeotype();
  void initChromosomesFromGenotypes();
  void initAllelesFromChromosomes();
  void initUsefulCross();
  void initDiffParents();
  void initTree();
  void initObligatoryLeaves();
  void initObjectiveGen();
  void initObjectiveCross();
  void initObjectivePop();
  void initBounds();
  bool isHomozygousBlock(int i, int p, int q) const;
  void printH() const;
  void printE() const;
  void printX() const;
  void printHX() const;
  void printHXE() const;
  void printInnerNodes() const;

  Genotype parseGenotype(int i) const;
  double parsePop(int i) const;
  double parseLogProb(int i) const;
  unsigned long parseGen(int i) const;
  void constructDAG();
  int getGap(int i) const;
  int getRelNode(int base_i, int backbone_i) const;
  int getAbsChromosome(int base_i, int backbone_k) const;
  int getAbsNode(int base_i, int backbone_i) const;
  int getRelChromosome(int base_k, int backbone_k) const;

protected:
  void printEdges(Node node, BoolNodeMap& visited, std::ostream& out) const;

public:
  IlpSolver(const Data* pData, const Options& options);
  ~IlpSolver();
  void init();
  SolverStatus solve(bool feasibility, int timeLimit);
  double getObjectiveValue() const;
};

inline int IlpSolver::getAbsChromosome(int i, int k) const
{
  if (k < 2*i)
    return k;
  else
    // KLOPT diT?
    return k + 2*(_options._bound - _options._fixedGen - i);
}

inline int IlpSolver::getAbsNode(int i, int j) const
{
  if (j < i)
    return j;
  else
    return j + _options._bound - _options._fixedGen - i;
}

inline int IlpSolver::getRelNode(int base_i, int backbone_i) const
{
  const int firstBackboneNode = _options._bound - _options._fixedGen;
  return base_i + backbone_i - firstBackboneNode;
}

inline int IlpSolver::getRelChromosome(int base_i, int backbone_k) const
{
  abort();
  // TODO
  const int firstBackboneChromosome = 2*(_options._bound - _options._fixedGen) - 1;
  return backbone_k - firstBackboneChromosome + 1 + 2*base_i;
}

inline int IlpSolver::getGap(int i) const
{
  assert(0 <= i && i < _options._bound);
  assert(_options._fixedGen > 1);

  return _options._bound - _options._fixedGen - i + 1;
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

inline void IlpSolver::printH() const
{
  const int m = _pData->getNumberOfLoci();
  for (int i = 0; i < _options._bound - 1; i++)
  {
    for (int p = 0; p < m-1; p++)
    {
      for (int q = p+1; q < m; q++)
      {
        std::cout << "// h[" << i << "][" << p << "][" << q << "] = "
          << _pCplex->getValue(_h[i][p][q-p-1]) << std::endl;
      }
    }
  }
}

inline void IlpSolver::printE() const
{
  const int m = _pData->getNumberOfLoci();
  for (int k = 0; k < 2*_options._bound - 1; k++)
  {
    for (int p = 0; p < m-1; p++)
    {
      for (int q = p+1; q < m; q++)
      {
        std::cout << "// e[" << k << "][" << p << "][" << q << "] = "
          << _pCplex->getValue(_e[k][p][q-p-1]) << std::endl;
      }
    }
  }
}

inline void IlpSolver::printX() const
{
  const int n = _pData->getParents().size();
  for (int k = 0; k < 2*_options._bound - 1; k++)
  {
    for (int i = 0; i < n+k/2; i++)
    {
      std::cout << "// x[" << k << "][" << i << "] = "
        << _pCplex->getValue(_x[k][i]) << std::endl;
    }
  }
}

inline void IlpSolver::printHX() const
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
          int val1 = _pCplex->getValue(_hx[k][i][p][q-p-1]);
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
            val2 = _pCplex->getValue(_h[i-n][p][q-p-1]) * _pCplex->getValue(_x[k][i]);
          }

          std::cout << "// hx[" << k << "][" << i << "]["
            << p << "][" << q << "] = " << val1 
            << (val1 == val2 ? "=" : "!=") << val2 << std::endl;
        }
      }
    }
  }
}

inline void IlpSolver::printHXE() const
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
          int val1 = _pCplex->getValue(_hxe[k][i][p][q-p-1]);
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
              val2 = _pCplex->getValue(_x[k][i]) * _pCplex->getValue(_e[k][p][q-p-1]);
            }
          }
          else
          {
            val2 = _pCplex->getValue(_h[i-n][p][q-p-1]) * 
              _pCplex->getValue(_x[k][i]) * _pCplex->getValue(_e[k][p][q-p-1]);
          }
          std::cout << "// hxe[" << k << "][" << i << "]["
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
    bool b0 = _pCplex->getValue(_p[2*i][l]) != 0;
    bool b1 = _pCplex->getValue(_p[2*i+1][l]) != 0;

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

  return _pCplex->getValue(_zb[i]);
}

inline unsigned long IlpSolver::parseGen(int i) const
{
  assert(0 <= i && i < _options._bound);
  return lround(_pCplex->getValue(_r[i]));
}

#endif /* ILPSOLVER_H_ */
