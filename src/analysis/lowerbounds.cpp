/*
 * lowerbounds.cpp
 *
 *  Created on: 3-may-2011
 *      Author: M. El-Kebir
 */

#include "lowerbounds.h"
#include <limits>
#include <ilcplex/ilocplex.h>
#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>

LowerBound::LowerBound(const Data* pData, bool verbose)
  : _pData(pData)
  , _verbose(verbose)
  , _popLB(std::numeric_limits<unsigned long>::max())
  , _crossLB(std::numeric_limits<unsigned long>::max())
  , _genLB(std::numeric_limits<unsigned long>::max())
{
  compute();
}

void LowerBound::compute()
{
  LinkageAnalysis popAnalysis(_pData);
  _popLB = 0;
  const LociPairList& unlinkedLoci = popAnalysis.getUnlinkedLoci(true);
	for (LociPairList::const_iterator it = unlinkedLoci.begin(); 
      it != unlinkedLoci.end(); it++)
	{
    if (it->first == it->second - 1)
    {
      _popLB += popAnalysis.getLinkagePopPair(true, *it).second;
    }
  }

  _crossLB = solveSetCover();
  _genLB = std::max(1., ceil(log(_crossLB)/log(2)));
  _popMaxLB = popAnalysis.getPopMax();

  if (_popLB < _popMaxLB)
    _popLB = _popMaxLB;
}

unsigned long LowerBound::solveSetCover()
{
  const GenotypeSet& parents = _pData->getParents();
  const int n = parents.size();
  const int m = _pData->getNumberOfLoci();

  IloEnv env;
  IloModel model(env);
  IloExpr obj(env);
  IloCplex* pCplex;

  IloBoolVarArray x(env, n);
  for (int i = 0; i < n; i++)
  {
    x[i] = IloBoolVar(env);
    obj += x[i];
  }

  IloExpr sum(env);
  for (int p = 0; p < m; p++)
  {
    int i=0;
    for (GenotypeSet::const_iterator it = parents.begin(); it != parents.end(); it++, i++)
    {
      const Genotype& C = *it;
      if (C(0,p) || C(1,p))
        sum += x[i];
    }
    model.add(sum >= 1);
    sum.clear();
  }
  sum.end();

  model.add(IloMinimize(env, obj));

  // UGLY hack, just to suppress the ILOG License Manager spamming
  int bakfd, newfd;
  fflush(stderr);
  bakfd = dup(2);
  newfd = open("/dev/null", O_WRONLY);
  dup2(newfd, 2);
  close(newfd);

  pCplex = new IloCplex(model);
  if (!_verbose)
  {
    pCplex->setOut(env.getNullStream());
    pCplex->setWarning(env.getNullStream());
    pCplex->setError(env.getNullStream());
  }
  else
  {
    pCplex->setOut(std::cerr);
    pCplex->setWarning(std::cerr);
    pCplex->setError(std::cerr);
  }

  fflush(stderr);
  dup2(bakfd, 2);
  close(bakfd);

  if (pCplex->solve())
  {
    _crossLB = pCplex->getObjValue();
    if (!_pData->isIdeotypeHomozygous()) _crossLB--;
  }


  delete pCplex;
  obj.end();
  model.end();

  return _crossLB;
}
