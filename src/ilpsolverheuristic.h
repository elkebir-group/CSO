/*
 * ilpsolverheuristic.h
 *
 *  Created on: 18-jun-2013
 *      Author: M. El-Kebir
 */

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>

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

  const int _n;
  const int _m;
  const int _nGen;
  const int _nCrs;
  BoolVarMatrix _x;

public:
  IlpSolverHeuristic(IloEnv env,
                     const int n,
                     const int m,
                     const int nGen,
                     const int nCrs,
                     BoolVarMatrix x);
  ~IlpSolverHeuristic() {}

protected:
  virtual void main();
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) IlpSolverHeuristic(*this));
  }
};

#endif // ILPSOLVERHEURISTIC_H
