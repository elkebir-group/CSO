/*
 * heuristicsolverfast.h
 *
 *  Created on: 09-apr-2009
 *      Author: M. El-Kebir
 */

#ifndef HEURISTICSOLVERFAST_H_
#define HEURISTICSOLVERFAST_H_

#include "solver.h"
#include "csodp.h"
#include "heuristicsolver.h"
#include "heuristiccorrectloci.h"
#include "heuristiccorrectalleles.h"
#include "heuristicindividualcost.h"
#include "heuristiclinkage.h"
#include "linkageanalysis.h"

class HeuristicSolverFast : public HeuristicSolver
{
private:
	int performStep(bool updateFlag);

public:
	HeuristicSolverFast(const Data* pData, DpTable* pTable, int maxIterations, int verbosity, 
		HeuristicType heuristicType, int nGenotypesToSelect, int nGenotypesLinkageAnalysis, bool selfingRestriction);
	~HeuristicSolverFast();
	virtual void solve();
};

#endif
