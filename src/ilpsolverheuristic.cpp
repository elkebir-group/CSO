#include "ilpsolverheuristic.h"

IlpSolverHeuristic::IlpSolverHeuristic(IloEnv env,
                                       const int n,
                                       const int m,
                                       const int nGen,
                                       const int nCrs,
                                       const BoolMatrix c,
                                       const IloExpr obj,
                                       BoolVarMatrix x,
                                       BoolVarMatrix p,
                                       IloNumVarArray allVar)

  : IloCplex::HeuristicCallbackI(env)
  , _n(n)
  , _m(m)
  , _nGen(nGen)
  , _nCrs(nCrs)
  , _c(c)
  , _obj(obj)
  , _x(x)
  , _p(p)
  , _allVar(allVar)
{
}

void IlpSolverHeuristic::main()
{
  // let's print x
  //std::cout << "gen = " << _nGen << std::endl;
  //std::cout << "crs = " << _nCrs << std::endl;

  //for (int k = 0; k < 2*_n; k++)
  //{
  //  for (int q = 0; q < _m; q++)
  //  {
  //    std::cout << "c[" << k << "][" << q << "] = "
  //              << _c[k][q] << std::endl;
  //  }
  //}

  //for (int k = 0; k < 2*_nCrs - 1; k++)
  //{
  //  int j = k / 2;
  //  //int offset = j >= _nCrs - _nGen ? 0 : _nGen - 2;
  //  //for (int i = 0; i < _n + j + offset; i++)
  //  //{
  //  //  std::cout << "x[" << k << "][" << i << "] = "
  //  //            << getValue(_x[k][i]) << std::endl;
  //  //}
  //  for (int q = 0; q < _m; q++)
  //  {
  //    std::cout << "p[" << k << "][" << q << "] = "
  //              << getValue(_p[k][q]) << std::endl;
  //  }
  //}
  //std::cout << std::endl;

  // let's construct something feasible:
  //
  // 1. first round x and p variables
  // 2. see whether this is feasible
  // 3. give back to cplex
  BoolVecMatrix chromosome;
  IntVector parent;
  for (int k = 0; k < 2*_nCrs - 1; k++)
  {
    // fix alleles
    chromosome.push_back(BoolVector());
    for (int q = 0; q < _m; q++)
    {
      chromosome.back().push_back(getValue(_p[k][q]) > 0.5);
    }

    // fix ancestor
    int j = k / 2;
    int offset = j >= _nCrs - _nGen ? 0 : _nGen - 2;
    int ancestor = -1;
    double ancestor_val = 0;
    for (int i = 0; i < _n + j + offset; i++)
    {
      double val = getValue(_x[k][i]);
      if (val > ancestor_val)
      {
        ancestor_val = val;
        ancestor = i;
      }
    }
    parent.push_back(ancestor);
  }

  if (isFeasibleSchedule(chromosome, parent))
  {
    //std::cout << "Found sth" << std::endl;
    IloNumVarArray new_var(getEnv());
    IloNumArray new_val(getEnv());

    for (int k = 0; k < 2*_nCrs - 1; k++)
    {
      // x vars
      int j = k / 2;
      int offset = j >= _nCrs - _nGen ? 0 : _nGen - 2;
      for (int i = 0; i < _n + j + offset; i++)
      {
        int bla = getFeasibility(_x[k][i]);
        int bla2 = CPX_IMPLIED_INTEGER_FEASIBLE;
        if (i == parent[k] && bla != bla2)
        {
          new_var.add(_x[k][i]);
          new_val.add(1);
        }
        if (i != parent[k] && bla != bla2)
        {
          new_var.add(_x[k][i]);
          new_val.add(0);
        }
      }

      // p vars
      for (int q = 0; q < _m; q++)
      {
        int bla = getFeasibility(_p[k][q]);
        int bla2 = CPX_IMPLIED_INTEGER_FEASIBLE;
        if (chromosome[k][q] && bla != bla2)
        {
          new_var.add(_p[k][q]);
          new_val.add(1);
        }
        if (!chromosome[k][q] && bla != bla2)
        {
          new_var.add(_p[k][q]);
          new_val.add(0);
        }
      }
    }

    //IloNumArray sol(getEnv(), _allVar.getSize());
    //getValues(sol, _allVar);
    //std::cout << "Number of variables: " << _allVar.getSize() << std::endl;
    //std::cout << "Number of variables: " << sol.getSize() << std::endl;
    //for (int i = 0; i < sol.getSize(); i++)
    //{
    //  const char* name = _allVar[i].getName();
    //  if (name)
    //  {
    //    std::cout << i << ": " << name << " = " << sol[i] << std::endl;
    //  }
    //  else
    //    std::cout << i << ": " << _allVar[i].getId() << " = " << sol[i] << std::endl;
    //}

    try
    {
      setBounds(new_var, new_val, new_val);
    }
    catch (IloCplex::Exception& e)
    {
      std::cout << e.getStatus() << std::endl;
    }

    // resolve
    double old_obj = getIncumbentObjValue();
    bool res = solve();
    if (res)
    {
      double new_obj = getValue(_obj);
      if (new_obj - old_obj > 1e-10)
      {
        IloNumArray sol(getEnv(), _allVar.getSize());
        getValues(sol, _allVar);
        for (int i = 0; i < sol.getSize(); i++)
        {
          if (fabs(sol[i]) < 1e-10) sol[i] = 0;
          //const char* name = _allVar[i].getName();
          //if (name)
          //{
          //  std::cout << i << ": " << name << " = " << sol[i] << std::endl;
          //}
          //else
          //  std::cout << i << ": " << _allVar[i].getId() << " = " << sol[i] << std::endl;
        }
        setSolution(_allVar, sol);

        std::cout << "New incumbent from " << old_obj << " to " << new_obj << std::endl;
      }
      else
      {
        std::cout << "Worse: " << new_obj << std::endl;
      }
    }
  }
}

void IlpSolverHeuristic::setSolutionLocal(const BoolVecMatrix& chromosome,
                                     const IntVector& parent)
{
  IloNumVarArray solutionVar(getEnv(), 0);
  IloNumArray solution;

  // 1. x
  add_x(parent, solutionVar, solution);
}

void IlpSolverHeuristic::add_x(const IntVector& parent,
                               IloNumVarArray& solutionVar,
                               IloNumArray& solution)
{
  IloNumArray solution_x;

  for (int k = 0; k < 2*_nCrs - 1; k++)
  {
    int j = k / 2;
    int offset = j >= _nCrs - _nGen ? 0 : _nGen - 2;
    for (int i = 0; i < _n + j + offset; i++)
    {
      solutionVar.add(_x[k][i]);

      if (i == parent[k]) solution_x.add(1);
      else solution_x.add(0);
    }
  }
  solution.add(solution_x);
}

void IlpSolverHeuristic::add_y(const BoolVecMatrix& chromosome,
                               const IntVector& parent,
                               IloNumVarArray& solutionVar,
                               IloNumArray& solution)
{
  IloNumArray solution_y;

  for (int k = 0; k < 2*_nCrs - 1; k++)
  {
    for (int q = 0; q < _m; q++)
    {
      //solutionVar.add(_y[k][q]);
    }
  }
}

bool IlpSolverHeuristic::isFeasibleSchedule(const BoolVecMatrix& chromosome,
                                            const IntVector& parent) const
{
  for (int k = 0; k < 2*_nCrs - 1; k++)
  {
    int anc = parent[k];
    if (anc < _n)
    {
      for (int q = 0; q < _m; q++)
      {
        bool allele = chromosome[k][q];
        bool anc_allele1 = _c[2 * anc][q];
        bool anc_allele2 = _c[2 * anc + 1][q];
        if (allele != anc_allele1 && allele != anc_allele2)
          return false;
      }
    }
    else
    {
      for (int q = 0; q < _m; q++)
      {
        bool allele = chromosome[k][q];
        bool anc_allele1 = chromosome[2 * (anc - _n)][q];
        bool anc_allele2 = chromosome[2 * (anc - _n) + 1][q];
        if (allele != anc_allele1 && allele != anc_allele2)
          return false;
      }
    }
  }

  return true;
}
