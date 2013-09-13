/*
 * ilpsolverheuristic.h
 *
 *  Created on: 18-jun-2013
 *      Author: M. El-Kebir
 */

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <vector>

#ifndef ILPSOLVERHEURISTIC_H
#define ILPSOLVERHEURISTIC_H

class IlpSolverHeuristic : public IloCplex::HeuristicCallbackI
{
private:
  typedef IloArray<IloBoolVarArray> BoolVarMatrix;
  typedef IloArray<BoolVarMatrix> BoolVar3Matrix;
  typedef IloArray<BoolVar3Matrix> BoolVar4Matrix;
  typedef IloArray<IloNumVarArray> NumVarMatrix;
  typedef IloArray<NumVarMatrix> NumVar3Matrix;
  typedef IloArray<IloIntVarArray> IntVarMatrix;
  typedef IloArray<IntVarMatrix> IntVar3Matrix;
  typedef IloArray<IloBoolArray> BoolMatrix;
  typedef std::vector<bool> BoolVector;
  typedef std::vector<BoolVector> BoolVecMatrix;
  typedef std::vector<int> IntVector;

  const int _n;
  const int _m;
  const int _nGen;
  const int _nCrs;
  const BoolMatrix _c;
  const IloExpr _obj;
  BoolVarMatrix _x;
  BoolVarMatrix _p;
  IloNumVarArray _allVar;

public:
  IlpSolverHeuristic(IloEnv env,
                     const int n,
                     const int m,
                     const int nGen,
                     const int nCrs,
                     const BoolMatrix c,
                     const IloExpr obj,
                     BoolVarMatrix x,
                     BoolVarMatrix p,
                     IloNumVarArray allVar);
  ~IlpSolverHeuristic() {}

protected:
  virtual void main();
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) IlpSolverHeuristic(*this));
  }

  bool isFeasibleSchedule(const BoolVecMatrix& chromosome,
                          const IntVector& parent) const;

  void setSolutionLocal(const BoolVecMatrix& chromosome,
                   const IntVector& parent);

  void add_x(const IntVector& parent,
             IloNumVarArray& solutionVar,
             IloNumArray& solution);

  void add_y(const BoolVecMatrix& chromosome,
             const IntVector& parent,
             IloNumVarArray& solutionVar,
             IloNumArray& solution);
};

#endif // ILPSOLVERHEURISTIC_H
