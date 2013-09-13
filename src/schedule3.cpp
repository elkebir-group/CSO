#include <iostream>
#include <sstream>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <string>
#include <boost/program_options.hpp>
#include <ilcplex/ilocplex.h>

#include "common/data.h"
#include "common/stopwatch.h"
#include "dp/dptable.h"
#include "dp/dptablehashmap.h"
#include "dp/dpitem.h"

using std::cerr;
using std::endl;
using std::cout;

namespace po = boost::program_options;

typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<BoolVarMatrix> BoolVar3Matrix;
typedef IloArray<IloNumVarArray> NumVarMatrix;

bool isTree(const DpItem& target)
{
  const DpItem* pParent1 = target.getParent1();
  const DpItem* pParent2 = target.getParent2();

  if (!pParent1 || !pParent2)
    return true;

  if (pParent1 == pParent2)
    return isTree(*pParent1);

  GenotypeSet ancestors1 = pParent1->getAncestors();
  GenotypeSet ancestors2 = pParent2->getAncestors();
  ancestors1.insert(*pParent1);
  ancestors2.insert(*pParent2);

  int res = 0;
  GenotypeSet::const_iterator it1 = ancestors1.begin(), it2 = ancestors2.begin();
  while (it1 != ancestors1.end() && it2 != ancestors2.end())
  {
    if (*it1 < *it2) ++it1;
    else if (*it2 < *it1) ++it2;
    else
    {
      res++;
      *it1++;
      *it2++;
    }
  }

  return res == 0;
}

bool isTarget(const int nLoci, const Genotype& genotype)
{
  // TODO: we need real ideotypes here
  const int target = (1 << nLoci) - 1;
  return genotype.getC0() == target || genotype.getC1() == target;
}

bool isNotFake(const int i,
               const bool tree,
               const int nParents,
               const IloCplex& cplex,
               const BoolVarMatrix& x)
{
  // check whether both chromosomes are from real parents
  int p1 = -1, p2 = -1;
  for (int j = 0; j < i + nParents; j++)
  {
    if (cplex.getValue(x[2*i][j]))
    {
      if (p1 != -1)
      {
        p1 = -1;
        break;
      }
      p1 = j;
    }
    if (cplex.getValue(x[2*i+1][j]))
    {
      if (p2 != -1)
      {
        p2 = -1;
        break;
      }
      p2 = j;
    }
  }

  //if (tree)
    return p1 != -1 || p2 != -1;
  //else
  //  return p1 != -1 && p2 != -1;
}

bool isNotFake2(const int i,
                const bool tree,
                const int nParents,
                const IloCplex& cplex,
                const BoolVar3Matrix& pp)
{
  // check whether both chromosomes are from real parents
  int p1 = -1, p2 = -1;
  for (int j = 0; j < i + nParents; j++)
  {
    if (cplex.getValue(pp[2*i][0][2*j]) || cplex.getValue(pp[2*i][0][2*j+1]) )
    {
      if (p1 != -1)
      {
        p1 = -1;
        break;
      }
      p1 = j;
    }
    if (cplex.getValue(pp[2*i+1][0][2*j]) || cplex.getValue(pp[2*i+1][0][2*j+1]) )
    {
      if (p2 != -1)
      {
        p2 = -1;
        break;
      }
      p2 = j;
    }
  }

  //if (tree)
    return p1 != -1 || p2 != -1;
  //else
  //  return p1 != -1 && p2 != -1;
}

const Genotype constructDAG(const Data* pData,
                            DpTable* pTable, 
                            const bool tree,
                            const bool multitarget,
                            const int bound,
                            const IloCplex& cplex,
                            const IloBoolVarArray& v, 
                            const BoolVarMatrix& x,
                            const BoolVarMatrix& p,
                            const BoolVarMatrix& d,
                            const IloIntVarArray& r)
{
  const int nLoci = pData->getNumberOfLoci();
  const GenotypeSet& parents = pData->getParents();
  const int nParents = parents.size();
  const DoubleMatrix& logRM = pData->getLogRM();

  GenotypePointerVector genotypeVector;
  GenotypeSet genotypeSet; 

  Genotype* pTarget = NULL;
  int targetGen = -1;

  // construct the parents
  for (GenotypeSet::const_iterator it = parents.begin();
    it != parents.end(); it++)
  {
    genotypeVector.push_back(new Genotype(*it));
    DpItem item(*it, NULL, NULL, 0, 0, 0, 0, 0, GenotypeSet(), false);
    pTable->updateItem(item);
  }

  for (int i = 0; i < bound; i++)
  {
    Genotype* pGenotype = NULL;

    // check whether a node is real
    if (cplex.getValue(v[i]) || isNotFake(i, tree, nParents, cplex, x))
    {
      // generate corresponding genotype
      int c0 = 0;
      int c1 = 0;

      double crossoverP1 = 0, crossoverP2 = 0;
      for (int l = 0; l < nLoci; l++)
      {
        bool b0 = cplex.getValue(p[2*i][l]) != 0;
        bool b1 = cplex.getValue(p[2*i+1][l]) != 0;
        c0 |= (int)b0 << (nLoci - l - 1);
        c1 |= (int)b1 << (nLoci - l - 1);

        if (l > 0 && cplex.getValue(d[2*i][l-1]))
          crossoverP1 += logRM[nLoci - l][nLoci - l - 1];
        if (l > 0 && cplex.getValue(d[2*i+1][l-1]))
          crossoverP2 += logRM[nLoci - l][nLoci - l - 1];
      }

      int gen = cplex.getValue(r[i]);
      pGenotype = new Genotype(c0, c1);

      if (genotypeSet.find(*pGenotype) == genotypeSet.end())
      {
        genotypeSet.insert(*pGenotype);

        DpItem item(*pGenotype, NULL, NULL, 
          crossoverP1, crossoverP2, gen, 0, 0, GenotypeSet(), false);
        pTable->updateItem(item);

        if (!pTarget && isTarget(nLoci, *pGenotype))
        {
          pTarget = pGenotype;
          targetGen = gen;
        }
      }
    }

    genotypeVector.push_back(pGenotype);
  }

  // make the edges
  for (int i = 0; i < bound; i++)
  {
    if (!cplex.getValue(v[i]) && genotypeVector[i+nParents] != pTarget)
      continue;

    int p1 = -1, p2 = -1;
    for (int j = 0; j < i + nParents; j++)
    {
      if (cplex.getValue(x[2*i][j]))
      {
        assert(p1 == -1);
        p1 = j;
      }
      if ((pTarget != genotypeVector[i+nParents] || multitarget) 
          && cplex.getValue(x[2*i+1][j]))
      {
        assert(p2 == -1);
        p2 = j;
      }
    }

    DpItem& item = pTable->getItem(*genotypeVector[i + nParents]);
    DpItem newItem = DpItem(item, 
      (p1 == -1 ? NULL : &pTable->getItem(*genotypeVector[p1])),
      (p2 == -1 ? NULL : &pTable->getItem(*genotypeVector[p2])), 
      item.getCrossoverP1(), 
      item.getCrossoverP2(), item.getGen(), 0, 0, 
      GenotypeSet(), false);

    if (newItem.getParent1() != &item && newItem.getParent2() != &item)
      item.update(newItem);
  }

  Genotype target = *pTarget;
  
  for (GenotypePointerVector::iterator it = genotypeVector.begin();
    it != genotypeVector.end(); it++)
  {
    delete *it;
  }

  pTable->getItem(target).updateAttributesFast(pTable, true);
  return target;
}

void printDAG(const Data* pData,
              const DpTable* pTable, 
              const DpItem& target,
              std::ostream& out)
{
  out << "digraph G {" << std::endl;
  target.printNodes(pTable, pData->getNumberOfLoci(), out);
  GenotypeSet genotypes;
  target.printEdges(out, pData->getNumberOfLoci(), genotypes);
  out << "}" << std::endl;
}

const DpItem* solveILP(const Data* pData,
                       DpTable* pTable,
                       bool verbose,
                       bool tree, 
                       bool diff_parents,
                       bool useful_cross,
                       bool multi_target)
{
  int _k = pData->getNumberOfLoci();

  const DoubleMatrix& logRM = pData->getLogRM();
  const DoubleMatrix& RM = pData->getRM();
  const GenotypeSet& parents = pData->getParents();
  int nof_parents = parents.size();

  const int nof_segs = 4;
  const double a[nof_segs+1] = {-7,-6,-5,-3,-0.2};

  //used only in multi target case - otherwise assumed single all 1's target
  const Genotype& ideotype = pData->getIdeotype();

  //int bound = _k; //(_k<nof_parents) ? _k : nof_parents;
  int bound = 5;

  IloEnv env;
  try
  {       
    //two targets - used only in multi-target case
    IloBoolArray t1(env,_k); 
    IloBoolArray t2(env,_k);
    for(int j=0; j<_k; j++)
    {
      t1[j] = (ideotype.getC0() >> _k-1-j) & 1;
      t2[j] = (ideotype.getC1() >> _k-1-j) & 1;
    }    

    IloArray<IloBoolArray> c(env, 2*nof_parents);

    // encode parent genotypes
    // c[i][j] = j-th bit of i-th chromosome
    int i=0;
    IloBoolArray homozygous(env,nof_parents);
    for (GenotypeSet::const_iterator it = parents.begin(); it != parents.end(); it++)
    {
      c[2*i] = IloBoolArray(env,_k);
      c[2*i+1] = IloBoolArray(env,_k);
      for(int j=0; j<_k; j++)
      {
        c[2*i][j] = (it->getC0() >> _k-1-j) & 1;
        c[2*i+1][j] = (it->getC1() >> _k-1-j) & 1;
      }
      homozygous[i]=it->isHomozygous();
      i++;
    }

    IloModel model(env);

    // x-Vars, encode edges
    // x[i][j] = chromosome i originates from node j
    BoolVarMatrix x(env, 2*bound);
    for(int i=0; i<bound; i++)
    {
      x[2*i] = IloBoolVarArray(env, nof_parents+i);
      x[2*i+1] = IloBoolVarArray(env, nof_parents+i);
      for(int j=0; j<nof_parents+i; j++)
      {
        x[2*i][j] = IloBoolVar(env);
        x[2*i+1][j] = IloBoolVar(env);
      }
    }

    BoolVar3Matrix input_bit(env, 2*bound); 
    for(int i=0; i<2*bound; i++)
    {           
      input_bit[i] = BoolVarMatrix(env,_k); 
      for(int j=0; j<_k; j++)
      {         
        input_bit[i][j] = IloBoolVarArray(env, 2*nof_parents + i - (i%2));
        for(int t=0; t<2*nof_parents+i-(i%2); t++) {
          std::stringstream ss;
          ss << "input_bit_";
          ss << i << "_" << j << "_" << t;        
          input_bit[i][j][t] = IloBoolVar(env,ss.str().c_str());
        }
      }      
    }

    // y-Vars, encodes bit origins
    // y[i][j] = { 0, if j-th bit of chromosome i originates 
    //                from the upper chromosome of the node defined by x
    //             1, if j-th bit of chromosome i originates
    //                from the lower chromosome of the node defined by x
    BoolVarMatrix y(env, 2*bound);
    for(int i=0; i<2*bound; i++)
    {
      y[i] = IloBoolVarArray(env, _k);
      for(int j=0; j<_k; j++)
        y[i][j] = IloBoolVar(env);
    }
  
    //TODO for trees: number of used leaves determines number of internal nodes

    //the bit variable
    BoolVarMatrix p(env, 2*bound);
    for(int i=0; i<2*bound; i++)
    {
      p[i] = IloBoolVarArray(env, _k);
      for(int t=0; t<_k; t++)
      {
        std::stringstream ss;
        ss << "p_" << i << "_" << t;
        p[i][t] = IloBoolVar(env,ss.str().c_str());
      }
    }

    // is locus different
    BoolVarMatrix bt(env, bound);
    for(int i=0; i<bound; i++)
    {
      bt[i] = IloBoolVarArray(env, _k);
      for(int t=0; t<_k; t++)
        bt[i][t] = IloBoolVar(env);
    }

    for(int i=0; i<bound; i++)
    {
      for(int j=0; j<_k; j++)
      {
        model.add(bt[i][j]>=p[2*i][j]-p[2*i+1][j]);  //TODO reicht das? keine <= Ungleichung notwendig?
        model.add(bt[i][j]>=p[2*i+1][j]-p[2*i][j]);
      }
    }

    // is node heterozygous
    IloBoolVarArray het(env,bound);
    for(int i=0; i<bound; i++)
    {
      het[i] = IloBoolVar(env);
      for(int j=0; j<_k; j++)
        model.add(het[i]>=bt[i][j]);
    }   

    //product var for heterozygous and x
    BoolVarMatrix hetx(env, 2*bound);
    for(int i=0; i<bound; i++)
    {
      hetx[2*i] = IloBoolVarArray(env, nof_parents+i);
      hetx[2*i+1] = IloBoolVarArray(env, nof_parents+i);
      for(int j=0; j<nof_parents; j++)
      {  
        hetx[2*i][j] = IloBoolVar(env);
        hetx[2*i+1][j] = IloBoolVar(env);
        if(homozygous[j]) {          
          model.add(hetx[2*i][j]==0);
          model.add(hetx[2*i+1][j]==0);
        }    
        else {      
          hetx[2*i][j] = x[2*i][j];
          hetx[2*i+1][j] = x[2*i+1][j];
        }
      }
      for(int j=nof_parents; j<nof_parents+i; j++)
      {
        hetx[2*i][j] = IloBoolVar(env);
        hetx[2*i+1][j] = IloBoolVar(env);
        model.add(hetx[2*i][j]<=het[j-nof_parents]);
        model.add(hetx[2*i][j]<=x[2*i][j]);
        model.add(hetx[2*i][j]>=het[j-nof_parents]+x[2*i][j]-1);

        model.add(hetx[2*i+1][j]<=het[j-nof_parents]);
        model.add(hetx[2*i+1][j]<=x[2*i+1][j]);
        model.add(hetx[2*i+1][j]>=het[j-nof_parents]+x[2*i+1][j]-1);      
      }
    }
  
    //TODO fuer fake nodes nix tun!

    //determine bit from parents
    IloExpr sumof1dep(env);
    for(int i=0; i<2*bound; i++)
    {
      for(int j=0; j<_k; j++)
      {
        for(int l=0; l<i-(i%2); l++)
        {
          model.add(input_bit[i][j][2*nof_parents+l]<=x[i][l/2 + nof_parents]);
          model.add(input_bit[i][j][2*nof_parents+l]<=p[l][j]);
          if((l%2) == 0) {
            model.add((input_bit[i][j][2*nof_parents+l]<=1-y[i][j]));
            model.add(input_bit[i][j][2*nof_parents+l]>=x[i][l/2+nof_parents]+p[l][j]+(1-y[i][j])-2);
          }
          else {
            model.add((input_bit[i][j][2*nof_parents+l]<=y[i][j]));               
            model.add(input_bit[i][j][2*nof_parents+l]>=x[i][l/2+nof_parents]+p[l][j]+y[i][j]-2);
          }
        }
        for(int l=0; l<2*nof_parents; l++)
        {
          model.add(input_bit[i][j][l]<=c[l][j]*x[i][l/2]);
          if((l%2) == 0) {
            model.add((input_bit[i][j][l]<=c[l][j]*(1-y[i][j])));
            model.add(input_bit[i][j][l]>=c[l][j]*(x[i][l/2]+(1-y[i][j])-1));
          }
          else {
            model.add((input_bit[i][j][l]<=c[l][j]*y[i][j]));
            model.add(input_bit[i][j][l]>=c[l][j]*(x[i][l/2]+y[i][j]-1));
          }
        }  
        for(int l=0; l<2*nof_parents+i-(i%2); l++)
          sumof1dep+=input_bit[i][j][l];    
        model.add(sumof1dep==p[i][j]);
        sumof1dep.clear();  
      }
    }
    sumof1dep.end();

    // //determine bit from parents
    // for(int i=0; i<2*_k; i++) {
    //   for(int i=0; i<2*bound; i++)
    //   {
    //     for(int t=0; t<_k; t++)
    //       {
    //         for(int j=0; j<i/2; j++)
    //     {
    //       model.add(p[i][t]>=1-y[i][t]+p[2*j][t]-1+x[i][nof_parents+j]-1);
    //       model.add(p[i][t]>=y[i][t]+p[2*j+1][t]-1+x[i][nof_parents+j]-1);
    //       model.add(p[i][t]<=p[2*j][t]+y[i][t]+1-x[i][nof_parents+j]);
    //       model.add(p[i][t]<=p[2*j+1][t]+1-y[i][t]+1-x[i][nof_parents+j]);
    //     }
    //         for(int j=0; j<nof_parents; j++)
    //     {
    //       // if(c[2*j][t]==1)
    //       //   model.add(p[i][t]>=-y[i][t]+x[i][j]);
    //       // else if(c[2*j][t]==0)
    //       //   model.add(p[i][t]<=y[i][t]+1-x[i][j]);

    //       // if(c[2*j+1][t]==1)
    //       //   model.add(p[i][t]>=y[i][t]-1+x[i][j]);
    //       // else if(c[2*j+1][t]==0)
    //       //   model.add(p[i][t]<=2-y[i][t]-x[i][j]);
              
    
    //       model.add(p[i][t]>=(1-y[i][t])*c[2*j][t]-1+x[i][j]);
    //       model.add(p[i][t]>=y[i][t]*c[2*j+1][t]-1+x[i][j]);
    //       //model.add(p[i][t]>=1-y[i][t]+c[2*j][t]-1+x[i][j]-1);
    //       //model.add(p[i][t]>=y[i][t]+c[2*j+1][t]-1+x[i][j]-1);
    //       model.add(p[i][t]<=c[2*j][t]+y[i][t]+1-x[i][j]);
    //       model.add(p[i][t]<=c[2*j+1][t]+1-y[i][t]+1-x[i][j]);    
    //     }
    //       }
    //   }
    // }
    //TODO in den c constraints vielleicht die Rollen von x und y tauschen?

    //BoolVarMatrix d(env, 2*_k);
    BoolVarMatrix d(env, 2*bound);
    for(int i=0; i<2*bound; i++)
    {
      d[i] = IloBoolVarArray(env, _k-1);
      for(int j=0; j<_k-1; j++)
        d[i][j] = IloBoolVar(env);
    }

    //did jth bit in string i result from crossing?
    //for(int i=0; i<2*_k; i++) {
    for(int i=0; i<2*bound; i++)
    {
      for(int j=1; j<_k; j++)
      {
        model.add(d[i][j-1]>=y[i][j]-y[i][j-1]);
        model.add(d[i][j-1]>=y[i][j-1]-y[i][j]);
      }
    }
   
    //product var of bt and x
    BoolVar3Matrix btx(env, 2*bound);
    for(int i=0; i<bound; i++)
    {
      btx[2*i] = BoolVarMatrix(env, nof_parents+i);
      btx[2*i+1] = BoolVarMatrix(env, nof_parents+i);
      for(int j=0; j<nof_parents; j++)
      {
        btx[2*i][j] = IloBoolVarArray(env,_k-1);
        btx[2*i+1][j] = IloBoolVarArray(env,_k-1);
        for(int t=0; t<_k-1; t++)
        {
          btx[2*i][j][t] = IloBoolVar(env);
          btx[2*i+1][j][t] = IloBoolVar(env);
          if(c[2*j][t+1]==c[2*j+1][t+1]) {
            model.add(btx[2*i][j][t]==0);
            model.add(btx[2*i+1][j][t]==0);  //TODO can these vars not be avoided?
          }
          else {
            model.add(btx[2*i][j][t]==x[2*i][j]); //TODO can these vars not be avoided?
            model.add(btx[2*i+1][j][t]==x[2*i+1][j]);
          }
        }
      }         
      for(int j=nof_parents; j<nof_parents+i; j++)
      {
        btx[2*i][j] = IloBoolVarArray(env,_k-1);
        btx[2*i+1][j] = IloBoolVarArray(env,_k-1);
        for(int t=0; t<_k-1; t++)
        {
          btx[2*i][j][t] = IloBoolVar(env);
          btx[2*i+1][j][t] = IloBoolVar(env);
          model.add(btx[2*i][j][t]<=x[2*i][j]);
          model.add(btx[2*i][j][t]<=bt[j-nof_parents][t+1]);
          model.add(btx[2*i][j][t]>=x[2*i][j]+bt[j-nof_parents][t+1]-1);

          model.add(btx[2*i+1][j][t]<=x[2*i+1][j]);
          model.add(btx[2*i+1][j][t]<=bt[j-nof_parents][t+1]);
          model.add(btx[2*i+1][j][t]>=x[2*i+1][j]+bt[j-nof_parents][t+1]-1);
        }
      }
    }
   
    //product var of bt and d and x (= product of btx and x)
    BoolVar3Matrix btxd(env, 2*bound);
    for(int i=0; i<bound; i++)
    {
      btxd[2*i] = BoolVarMatrix(env, nof_parents+i);
      btxd[2*i+1] = BoolVarMatrix(env, nof_parents+i);
      for(int j=0; j<nof_parents; j++)
      {
        btxd[2*i][j] = IloBoolVarArray(env,_k-1);
        btxd[2*i+1][j] = IloBoolVarArray(env,_k-1);
        for(int t=0; t<_k-1; t++)
        {
          btxd[2*i][j][t] = IloBoolVar(env);
          btxd[2*i+1][j][t] = IloBoolVar(env);
          if(c[2*j][t+1]==c[2*j+1][t+1]) {
            model.add(btxd[2*i][j][t]==0);
            model.add(btxd[2*i+1][j][t]==0); //TODO can these vars not be avoided?
          }
          else {
            model.add(btxd[2*i][j][t]<=x[2*i][j]);
            model.add(btxd[2*i][j][t]<=d[2*i][t]);
            model.add(btxd[2*i][j][t]>=x[2*i][j]+d[2*i][t]-1);

            model.add(btxd[2*i+1][j][t]<=x[2*i+1][j]);
            model.add(btxd[2*i+1][j][t]<=d[2*i+1][t]);
            model.add(btxd[2*i+1][j][t]>=x[2*i+1][j]+d[2*i+1][t]-1);
          }
        }
      }
      for(int j=nof_parents; j<nof_parents+i; j++)
      {
        btxd[2*i][j] = IloBoolVarArray(env,_k-1);
        btxd[2*i+1][j] = IloBoolVarArray(env,_k-1);
        for(int t=0; t<_k-1; t++)
        {
          btxd[2*i][j][t] = IloBoolVar(env);
          btxd[2*i+1][j][t] = IloBoolVar(env);      
          model.add(btxd[2*i][j][t]<=btx[2*i][j][t]); //TODO or better directly in terms of bt and x?
          model.add(btxd[2*i][j][t]<=d[2*i][t]);
          model.add(btxd[2*i][j][t]>=btx[2*i][j][t]+d[2*i][t]-1);
          
          model.add(btxd[2*i+1][j][t]<=btx[2*i+1][j][t]);
          model.add(btxd[2*i+1][j][t]<=d[2*i+1][t]);
          model.add(btxd[2*i+1][j][t]>=btx[2*i+1][j][t]+d[2*i+1][t]-1);
        }
      }
    }      
   
    //separable z var
    IloNumVarArray zb(env, bound, -IloInfinity, 0.0, IloNumVar::Float);
    IloExpr sumoflogprobs(env);
    for(int i=0; i<bound; i++) 
    {      
      for(int j=0; j<nof_parents+i; j++) {

        // MEK: it chrashes here if parent is homozygous, apparently 0.5 * 0 does not work
        if (!(j < nof_parents && homozygous[j]))
        {
          sumoflogprobs+=log(0.5)*hetx[2*i][j];
          sumoflogprobs+=log(0.5)*hetx[2*i+1][j];
        }
        for(int t=0; t<_k-1; t++)
        {           
          double pt = RM[_k-t-1][_k-t-2];
          sumoflogprobs+=log(1-pt)*btx[2*i][j][t];
          sumoflogprobs+=log(1-pt)*btx[2*i+1][j][t];
          sumoflogprobs+=log(pt/(1.-pt))*btxd[2*i][j][t];
          sumoflogprobs+=log(pt/(1.-pt))*btxd[2*i+1][j][t];
        }       
      }
      model.add(zb[i]==sumoflogprobs);
      sumoflogprobs.clear();
    }
    sumoflogprobs.end();
    
    NumVarMatrix lambda(env, bound);
    IloExpr sumoflambdas(env);
    for(int i=0; i<bound; i++)
    {
      lambda[i] = IloNumVarArray(env, nof_segs+1, 0., 1., IloNumVar::Float);
           
      for(int j=0; j<nof_segs+1; j++)
        sumoflambdas+=lambda[i][j];
      model.add(sumoflambdas==1);
      sumoflambdas.clear();

      for(int j=0; j<nof_segs+1; j++)
        sumoflambdas+=a[j]*lambda[i][j];
      model.add(zb[i]==sumoflambdas);
      sumoflambdas.clear();    
    }
    sumoflambdas.end();     

    //either this block....************************************

    // //IloBoolVarArray v(env,_k);
    // IloBoolVarArray v(env,bound);
    // //count real internal nodes
    // for(int i=0; i<bound; i++)
    // {
    //   v[i] = IloBoolVar(env);
    //   for(int j=0; j<_k; j++) {
    //   if(multi_target)
    //     {
    //       if(t1[j])
    //         model.add(v[i]>=1-p[2*i][j]);
    //       else
    //         model.add(v[i]>=p[2*i][j]);
      
    //       if(t2[j])
    //         model.add(v[i]>=1-p[2*i+1][j]);
    //       else
    //         model.add(v[i]>=p[2*i+1][j]);

    //     } else
    //     {
    //       model.add(v[i]>=1-p[2*i][j]);
    //     }
    //   }
    // }
    
    // //only fake nodes after a fake node  -- seems to make it worse
    // // for(int i=1; i<bound; i++)       
    // //   model.add(v[i]<=v[i-1]);        

    // //each chromosome (except for last) has exactly one parent node
    // IloExpr sumofinedges(env);
    // //for(int i=0; i<2*_k; i++) {
    // for(int i=0; i<2*bound-1; i++) //for bound-1 much slower for tree (badexmpl) with no bounds
    //   {                            //but much faster if bound k=4 used
    //   for(int j=0; j<nof_parents; j++)
    //     sumofinedges+=x[i][j];
    //   for(int j=0; j<i/2; j++)
    //     sumofinedges+=x[i][j+nof_parents];
    //   if(i%2==0)
    //   model.add(sumofinedges==1);
    //   if(i%2==1) {
    //   if(multi_target)
    //     model.add(sumofinedges==1);
    //   else
    //     model.add(sumofinedges==v[i/2]); //TODO why slower?
    //   }
    //   sumofinedges.clear();
    // }
    // sumofinedges.end();

    // //feasible -- needed ? since minimizing...
    // for(int j=0; j<_k; j++) {
    //   if(multi_target)
    //   {
    //     model.add(p[2*bound-2][j]==t1[j]);
    //     model.add(p[2*bound-1][j]==t2[j]);
    //   } else
    //   {
    //     model.add(p[2*bound-2][j]==1);
    //   }      
    // }
    // //model.add(v[bound-1]==0); seems to make it slower


    //or this block.. ********************************************

    //real node cannot succeed fake node
    IloBoolVarArray v(env,bound);
    for(int i=1; i<bound; i++)       
      model.add(v[i]<=v[i-1]);    

    IloExpr sumofinedges1(env);
    IloExpr sumofinedges2(env);
    for(int i=0; i<bound; i++)
    {
      for(int j=0; j<nof_parents; j++) {
        sumofinedges1+=x[2*i][j];
        sumofinedges2+=x[2*i+1][j]; 
      }
      for(int j=0; j<i; j++) {
        sumofinedges1+=x[2*i][j+nof_parents];
        sumofinedges2+=x[2*i+1][j+nof_parents];
      }
      model.add(sumofinedges1<=1);
      //model.add(sumofinedges2<=1);
      model.add(sumofinedges1==sumofinedges2);
      model.add(v[i]==sumofinedges1);
      sumofinedges1.clear();
      sumofinedges2.clear();
    }
    sumofinedges1.end();
    sumofinedges2.end();

    //node is target
    IloBoolVarArray t(env,bound);
    for(int i=0; i<bound; i++)
    {
      t[i] = IloBoolVar(env);

      //fake node cannot be target
      model.add(t[i]<=v[i]);

      for(int j=0; j<_k; j++) 
      {         
        if(multi_target) {
          if(t1[j])
            model.add(t[i]<=p[2*i][j]);
          else
            model.add(t[i]<=1-p[2*i][j]);
      
          if(t2[j])
            model.add(t[i]<=p[2*i+1][j]);
          else
            model.add(t[i]<=1-p[2*i+1][j]);
        } else {    
          model.add(t[i]<=p[2*i][j]);
        }
      }
    }
    IloExpr sumoftargets(env);
    for(int i=0; i<bound; i++)   
      sumoftargets+=t[i];
    model.add(sumoftargets==1); //TODO or <= ?
    sumoftargets.end();
      

    // *****************************************************************

    if(useful_cross)
    {
      //collect at least one '1' between two crossings in opt
      IloExpr sumofcollected(env);
      for(int i=0; i<bound; i++)
      {
        for(int j=0; j<_k-2; j++)
          for(int l=j+1; l<_k-1; l++)
          {
            for(int t=j+1; t<l+1; t++)
              sumofcollected+=p[i][t];
            model.add(d[i][j]+d[i][l]-sumofcollected<=1);
            sumofcollected.clear();
          }
        for(int j=0; j<_k-1; j++) {
          for(int t=0; t<j+1; t++)
            sumofcollected+=p[i][t];
          model.add(d[i][j]-sumofcollected<=0);
          sumofcollected.clear();
        }      
        for(int j=0; j<_k-1; j++) {
          for(int t=j+1; t<_k; t++)
            sumofcollected+=p[i][t];
          model.add(d[i][j]-sumofcollected<=0);
          sumofcollected.clear();
        }
      }
      sumofcollected.end();
    }
    
    if(diff_parents)
    {
      for(int i=0; i<bound; i++)
      {
        for(int j=0; j<nof_parents; j++)
          model.add(x[2*i][j]+x[2*i+1][j]<=1);
        for(int j=0; j<i; j++)
          model.add(x[2*i][j+nof_parents]+x[2*i+1][j+nof_parents]<=1);
      }
    }
   
    if(tree) 
    {
      //every real node is used at most once as parent (=tree)
      IloExpr sumofoutedges(env);
      for(int i=0; i<nof_parents; i++)
      {
        for(int j=0; j<2*bound; j++)
        {
          sumofoutedges+=x[j][i];
        }
        model.add(sumofoutedges<=1);
        sumofoutedges.clear();
      }
      for(int i=0; i<bound-1; i++)
      {
        for(int j=i+1; j<bound; j++)
        {
          sumofoutedges+=x[2*j][nof_parents+i];
          sumofoutedges+=x[2*j+1][nof_parents+i];
        }
        model.add(sumofoutedges<=1);
        sumofoutedges.clear();
      }
      sumofoutedges.end();
    }

    //outdegree of obligatory leaves at least one
    bool* obligatory = new bool[nof_parents];
    for(int i=0; i<nof_parents; i++)
      obligatory[i]=false;
   
    for(int j=0; j<_k; j++)
    {
      int i=0;
      bool unique_covered = false;
      int obl_parent = -1;
      for (GenotypeSet::const_iterator it = parents.begin(); it != parents.end(); it++)
      {
        if(c[2*i][j] || c[2*i+1][j])
        {
          if(unique_covered) {
            unique_covered=false;
            break;
          }
          else {
            unique_covered=true;
            obl_parent = i;
          }
        }
        i++;
      }
      if(unique_covered)
      {     
        obligatory[obl_parent]=true;
      }
    }

    IloExpr sumofoutedges(env);
    for(int i=0; i<nof_parents; i++)
    {
      if(!obligatory[i])
        continue;
      for(int j=0; j<2*bound; j++)
        {    
          sumofoutedges+=x[j][i];
        }
      model.add(sumofoutedges>=1);       
      sumofoutedges.clear();
    }    
    sumofoutedges.end();
    delete obligatory;

    IloIntVarArray r(env, bound, 0, bound);
    for(int i=0; i<bound; i++)
    {
      for(int j=0; j<i; j++)
      {
        model.add(r[i]>=r[j]+1-(bound*(1-x[2*i][nof_parents+j]))-(bound*(1-v[j])));
        model.add(r[i]>=r[j]+1-(bound*(1-x[2*i+1][nof_parents+j]))-(bound*(1-v[j])));
      }
      for(int j=0; j<nof_parents; j++)
      {
        //TODO not needed! >=1 sufficient!!!
        model.add(r[i]>=1-(bound*(1-x[2*i][j])));
        model.add(r[i]>=1-(bound*(1-x[2*i+1][j])));
      }
    }
   
    IloIntVar depth(env, 0, bound);
    for(int i=0; i<bound; i++)
      model.add(depth>=r[i]);  

    IloExpr obj(env);
  
    if (pData->getCostNode() != 0)
      for(int j=0; j<bound; j++)
        obj += pData->getCostNode() * v[j]; // v[j] = true iff j is a true node

    // if (pData->getCostCrossover() != 0) {     
    //   IloExpr sumofcross(env);
    //   for(int j=1; j<_k; j++)
    //       {
    //         for(int i=0; i<2*bound; i++)
    //           {
    //             sumofcross += logRM[_k - j][_k - j - 1] * d[i][j-1]; //TODO remove fake chromosomes from some? 
    //           }
    //         obj -= pData->getCostCrossover() * sumofcross;  
    //         sumofcross.clear();
    //       } 
    //   sumofcross.end();      
    // }

    IloExpr sumofpopsize(env);
    for(int i=0; i<bound; i++)
    {
      for(int j=0; j<nof_segs+1; j++)
    //sumofpopsize+=(log(1-pData->getGamma())/log(1-exp(a[j])))*lambda[i][j]; 
    sumofpopsize+=(log(1-pData->getGamma())/log(1-exp(a[j])))*lambda[i][j];  //TODO where is gamma?
    }
    obj += pData->getCostCrossover() * sumofpopsize;
    sumofpopsize.clear();
    sumofpopsize.end();

    if (pData->getCostGen() != 0)
      obj += pData->getCostGen() * depth;

    //IloIntVar zz(env);
    //model.add(zz==obj);
    //model.add(zz<=16);

    model.add(IloMinimize(env,obj));

    // UGLY hack, just to suppress the ILOG License Manager spamming
    int bakfd, newfd;
    fflush(stderr);
    bakfd = dup(2);
    newfd = open("/dev/null", O_WRONLY);
    dup2(newfd, 2);
    close(newfd);

    IloCplex cplex(model);

    fflush(stderr);
    dup2(bakfd, 2);
    close(bakfd);

    if (!verbose)
    {
      cplex.setOut(env.getNullStream());
      cplex.setWarning(env.getNullStream());
      cplex.setError(env.getNullStream());
    }
    else
    {
      cplex.setOut(std::cerr);
      cplex.setWarning(std::cerr);
      cplex.setError(std::cerr);
    }

    //cplex.exportModel("letssee.lp");
    if( !cplex.solve() )
    {
      env.error() << "Failed to optimize LP" << endl;
      return NULL;
    }
    cout << cplex.getObjValue() << endl;
    //exit(0);

    Genotype targetGenotype = constructDAG(pData, pTable, tree, multi_target, bound, cplex, v, x, p, d, r);
    return &pTable->getItem(targetGenotype);
  }
  catch (IloException& ex)
  {
    cerr << "Error: " << ex << endl;
    return NULL;
  }
}


int main(int argc, char* argv[])
{
  bool verbose = false;
  bool tree = false;
  std::string inputFileName;
  int formulation=0;
  bool diff_parents = false;
  bool useful_cross = false;
  bool multi_target = false;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Help message")
    ("verbose,v", po::value<bool>(&verbose)->zero_tokens(), "Verbose CPLEX output")
    ("tree,t", po::value<bool>(&tree)->zero_tokens(), "Enforce trees")
    ("input-file,i", po::value<std::string>(&inputFileName), "Input file name")
    ("formulation,f", po::value<int>(&formulation), "Formulation")
    ("different-parents,p", po::value<bool>(&diff_parents)->zero_tokens(), "Cuts: Different parents") //note: does not make sense for two identical targets
    ("useful-crossovers,c", po::value<bool>(&useful_cross)->zero_tokens(), "Cuts: Useful crossovers only")
    ("multi-target,m", po::value<bool>(&multi_target)->zero_tokens(), "Multi target");

  try
  {
    po::positional_options_description p;
    p.add("input-file", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
      std::cout << "Usage: " << argv[0] << " [options]\n" << desc << std::endl;
      return 0;
    }

    if (!vm.count("input-file"))
    {
      std::cerr << "Missing input file" << std::endl;
      return 1;
    }
  }
  catch (std::exception&)
  {
    std::cerr << "Usage: " << argv[0] << " [options]\n" << desc << std::endl;
    return 1;
  }

  Data* pData = Data::create(inputFileName.c_str());
  if (!pData) return 1;

  if (verbose)
  {
    std::cout << "// RM" << std::endl;
    pData->printRM();
    std::cout << "// log RM" << std::endl;
    pData->printLogRM();
  }

  StopWatch stopWatch;
  DpTableHashMap table(pData->getNumberOfLoci());

  stopWatch.start();
  const DpItem* pTarget = solveILP(pData, &table, verbose, tree, diff_parents, useful_cross, multi_target);
  stopWatch.stop();

  if (pTarget)
  {
    printDAG(pData, &table, *pTarget, std::cout);

    fprintf(stderr, "\"%s\",%.3f,%lu,%lu,%.2f,%.2f,\"%s\"\n", 
      inputFileName.c_str(), stopWatch.getElapsedTime(),
      pTarget->getGen(), pTarget->getCumCross(), pTarget->getCumCrossover(),
      pData->getCost(pTarget->getCumCrossover(), pTarget->getGen(), pTarget->getCumCross()),
      isTree(*pTarget) ? "tree" : "DAG");

    delete pData;
    return 0;
  }
  else
  {
    fprintf(stderr, "\"%s\",%.3f,-,-,-,-,-\n", inputFileName.c_str(), stopWatch.getElapsedTime());

    delete pData;
    return 1;
  }
}
