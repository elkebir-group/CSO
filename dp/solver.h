/*
 * solver.h
 *
 *  Created on: 18-feb-2009
 *      Author: M. El-Kebir
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "csodp.h"
#include "data.h"
#include "genotype.h"
#include "dpitem.h"
#include <math.h>
#include "dptable.h"
#include <string>

class Solver
{
private:
	int performStep(bool updateFlag);

protected:
	const Data* _pData;
	DpTable* _M;
	GenotypeSet _currentGenotypes;
	GenotypeSet _oldGenotypes;
	const int _maxIterations;
	double _ideotypeCost;
	int _iteration;
	const int _verbosity;
	std::string _resDAG;
	
	int cleanUp(GenotypeSet& genotypeSet);
	int cleanUp(DpItemSet& dpItemSet);
	virtual void cross(const Genotype& genotype1, const Genotype& genotype2, 
		DpItemSet& newItemSet, bool updateFlag);
	void updateSet(const DpItem& newItem, const double newItemScore, DpItemSet& itemSet);

public:
	static Solver* _DP;

	Solver(const Data* pData, DpTable* pTable, int maxIterations, int verbosity);
	virtual ~Solver();
	virtual void solve();
	double getScore(const DpItem& item) const;
	void printDAG(const DpItem& item, std::ostream& out) const;
	const std::string& getResDAG() const;
	virtual std::string getMethodName() const;
};

inline int Solver::cleanUp(GenotypeSet& genotypeSet)
{
	int i = 0;

	GenotypeSet::iterator it = genotypeSet.begin();
	GenotypeSet::iterator newIt = it;
	while (it != genotypeSet.end())
	{
		newIt++;

		if (getScore(_M->getItem(*it)) > _ideotypeCost)
		{
			i++;
			genotypeSet.erase(it);
		}

		it = newIt;
	}

	return i;
}

inline int Solver::cleanUp(DpItemSet& dpItemSet)
{
	int i = 0;
	
	DpItemSet::iterator it = dpItemSet.begin();
	DpItemSet::iterator newIt = it;
	while (it != dpItemSet.end())
	{
		newIt++;

		if (getScore(*it) > _ideotypeCost)
		{
			i++;
			dpItemSet.erase(it);
		}

		it = newIt;
	}
	
	return i;
}

inline const std::string& Solver::getResDAG() const
{
	return _resDAG;
}

#endif /* SOLVER_H_ */
