#include "ilpsolverheuristic.h"

IlpSolverHeuristic::IlpSolverHeuristic(IloEnv env,
                                       const int n,
                                       const int m,
                                       const int nGen,
                                       const int nCrs,
                                       BoolVarMatrix x)
  : IloCplex::HeuristicCallbackI(env)
  , _n(n)
  , _m(m)
  , _nGen(nGen)
  , _nCrs(nCrs)
  , _x(x)
{
}

void IlpSolverHeuristic::main()
{
  // TODO
  if (hasIncumbent())
  {
    //std::cout << "Hoi: " << getIncumbentObjValue() << std::endl;
  }
  else
  {
    // let's print x
    for (int k = 0; k < 2*_nCrs - 1; k++)
    {
      int j = k / 2;
      int offset = j >= _nCrs - _nGen ? 0 : _nGen - 2;
      for (int i = 0; i < _n + j + offset; i++)
      {
        std::cout << "x[" << k << "][" << i << "] = "
                  << getValue(_x[k][i]) << std::endl;
      }
    }
    std::cout << std::endl;
  }
}
