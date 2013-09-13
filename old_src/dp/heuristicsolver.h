/*
 * heuristicsolver.h
 *
 *  Created on: 08-apr-2009
 *      Author: M. El-Kebir
 */

#ifndef HEURISTICSOLVER_H_
#define HEURISTICSOLVER_H_

#include "solver.h"
#include "csodp.h"
#include "../analysis/linkageanalysis.h"

class HeuristicSolver : public Solver
{
private:
	int performStep(bool updateFlag);
	void applySortingHeuristic(DpItemVector& newItemVector, int target);

protected:
	HeuristicType _heuristicType;
	int _nGenotypesToSelect;
	LinkageAnalysis _linkageAnalysis;
	int _nGenotypesLinkageAnalysis;
	bool _selfingRestriction;

	DpItemSet applyHeuristic(const DpItemSet& newItemSet, int target);
	void applyLinkageAnalysisUnlinked(const DpItemVector& newItemVector, const LociPairList& lociList, 
		const int targetChromosome, DpItemSet& itemsToKeep) const;
	void applyLinkageAnalysisWeaklyLinked(const DpItemVector& newItemVector, const LociPairList& lociList, 
		const int targetChromosome, DpItemSet& itemsToKeep) const;
	int applyInitialParentSelfing();
	virtual void cross(const Genotype& genotype1, const Genotype& genotype2, 
		DpItemSet& newItemSet, bool updateFlag);

public:
	HeuristicSolver(const Data* pData, DpTable* pTable, int maxIterations, int verbosity, 
		HeuristicType heuristicType, int nGenotypesToSelect, int nGenotypesLinkageAnalysis,
		bool selfingRestriction);
	virtual ~HeuristicSolver();
	virtual void solve();
	virtual std::string getMethodName() const;
};

#endif
