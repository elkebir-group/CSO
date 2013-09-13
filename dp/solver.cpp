/*
 * solver.cpp
 *
 *  Created on: 18-feb-2009
 *      Author: M. El-Kebir
 */

#include "solver.h"
#include "heuristicindividualcost.h"
#include "heuristiccorrectloci.h"
#include "heuristiccorrectalleles.h"
#include "heuristiclinkage.h"
#include <iostream>
#include <algorithm>
#include <sstream>

Solver* Solver::_DP = NULL;

Solver::Solver(const Data* pData, DpTable* pTable, int maxIterations, int verbosity)
	: _pData(pData)
	, _M(pTable)
	, _currentGenotypes(pData->getParents())
	, _oldGenotypes()
	, _maxIterations(maxIterations)
	, _ideotypeCost(DBL_MAX)
	, _iteration(0)
	, _verbosity(verbosity)
{
	const GenotypeSet& parents = pData->getParents();

	for (GenotypeSet::const_iterator it = parents.begin();
		it != parents.end();
		it++)
	{
		_M->updateItem(DpItem(*it, false));
	}

	_DP = this;
}

Solver::~Solver()
{
}

void Solver::solve()
{
	bool updateFlag = true;
	int newIndividuals = (int) _currentGenotypes.size();
	for (int _iteration = 0; _iteration < _maxIterations && newIndividuals; _iteration++)
	{
		if (_verbosity > 0)
		{
			std::cout << "\n// Iteration: " << _iteration + 1 << "\t Number of (current/old) individuals: " 
				<< (int) _currentGenotypes.size() << '/' << (int) _oldGenotypes.size() << std::endl;
		}
		
		updateFlag = !updateFlag;
		newIndividuals = performStep(updateFlag);
	}

	if (_verbosity > 0)
	{
		std::cout << "\n// Size of _oldGenotypes " << (int) _oldGenotypes.size() << " genotypes." << std::endl;
		std::cout << "// Number of generated genotypes: " << _M->getUseCount() << std::endl;
	}
}

int Solver::performStep(bool updateFlag)
{
	/*for (GenotypeSet::iterator it = _oldGenotypes.begin(); it != _oldGenotypes.end(); it++)
	{
		DpItem& item = _M->getItem(*it);
		if (item.getUpdateFlag() != updateFlag)
		{
			std::cout << "ERROR - invalid updateflag\n";
		}
		if ((item.getParent1() && item.getCumPop() < item.getParent1()->getCumPop()) || 
			(item.getParent2() && item.getCumPop() < item.getParent2()->getCumPop()))
		{
			std::cout << "ERROR - invalid popSize\t";
			item.printItem(Data::getInstance()->getNumberOfLoci());
		}
	}*/

	DpItemSet newItemSet;

	if (_verbosity > 0)
		std::cout << "// Performing " << _currentGenotypes.size() * (_currentGenotypes.size() + 1) / 2 
			<< " crossings among current individuals." << std::endl;

	/* perform crossings among current individuals */
	if (_pData->getAllowSelfing())
	{
		for (GenotypeSet::iterator it1 = _currentGenotypes.begin(); it1 != _currentGenotypes.end(); it1++)
		{
			for (GenotypeSet::iterator it2 = it1; it2 != _currentGenotypes.end(); it2++)
			{
				cross(*it1, *it2, newItemSet, updateFlag); // O(4^m)
			}
		}
	}
	else
	{
		for (GenotypeSet::iterator it1 = _currentGenotypes.begin(); it1 != _currentGenotypes.end(); it1++)
		{
			GenotypeSet::iterator it2 = it1;
			for (it2++; it2 != _currentGenotypes.end(); it2++)
			{
				cross(*it1, *it2, newItemSet, updateFlag); // O(4^m)
			}
		}
	}

	int erasedCount;
	if ((erasedCount = cleanUp(_oldGenotypes)))
	{
		if (_verbosity > 0)
			std::cout << "// Erased " << erasedCount << " old sub-optimal genotypes." << std::endl;
	}
	if ((erasedCount = cleanUp(_currentGenotypes)))
	{
		if (_verbosity > 0)
			std::cout << "// Erased " << erasedCount << " current sub-optimal genotypes." << std::endl;
	}

	if (_verbosity > 0)
		std::cout << "// Performing " << _currentGenotypes.size() * _oldGenotypes.size() 
			<< " crossings among current and old individuals." << std::endl;

	/* perform crossings among current and old individuals */
	for (GenotypeSet::iterator it1 = _currentGenotypes.begin(); it1 != _currentGenotypes.end(); it1++)
	{
		for (GenotypeSet::iterator it2 = _oldGenotypes.begin(); it2 != _oldGenotypes.end(); it2++)
		{
			cross(*it1, *it2, newItemSet, updateFlag); // O(2^m)
		}
	}

	int cheaper = 0;
	for (DpItemSet::iterator it = newItemSet.begin(); it != newItemSet.end(); it++)
	{
		cheaper += _currentGenotypes.erase(*it);
		cheaper += _oldGenotypes.erase(*it);
		
		_M->updateItem(*it);
	}

	if (_verbosity > 0)
	{
		std::cout << "// Number of new individuals: " << newItemSet.size() - cheaper << std::endl;
		std::cout << "// Number of cheaper individuals: " << cheaper << std::endl;
	}

	set_union(_oldGenotypes.begin(), _oldGenotypes.end(),
		_currentGenotypes.begin(), _currentGenotypes.end(),
		inserter(_oldGenotypes, _oldGenotypes.begin()));

	if ((erasedCount = cleanUp(_oldGenotypes)))
	{
		if (_verbosity > 0)
			std::cout << "// Erased " << erasedCount << " old sub-optimal genotypes." << std::endl;
	}

	_currentGenotypes.clear();
	_currentGenotypes.insert(newItemSet.begin(), newItemSet.end());

	for (GenotypeSet::reverse_iterator it = _currentGenotypes.rbegin(); it != _currentGenotypes.rend(); it++)
	{
		_M->getItem(*it).updateAttributesFast(_M, !updateFlag);
	}

	for (GenotypeSet::reverse_iterator it = _oldGenotypes.rbegin(); it != _oldGenotypes.rend(); it++)
	{
		_M->getItem(*it).updateAttributesFast(_M, !updateFlag);
	}

	if ((erasedCount = cleanUp(_currentGenotypes)))
	{
		if (_verbosity > 0)
			std::cout << "// Erased " << erasedCount << " current sub-optimal genotypes." << std::endl;
	}

	return (int) newItemSet.size();
}

void Solver::updateSet(const DpItem& newItem, const double newItemScore, DpItemSet& itemSet)
{
	DpItemSet::iterator itemIt = itemSet.find(newItem);
	if (itemIt == itemSet.end())
	{
		itemSet.insert(newItem);
	}
	else if (newItemScore < getScore(*itemIt))
	{
		// bug in g++
		DpItem& test = const_cast<DpItem&>(*itemIt);
		test.update(newItem);
	}

	if (newItem == _pData->getIdeotype())
	{
		if (_verbosity > 1)
		{
			//newItem.printItem(_pData->getNumberOfLoci());
			printDAG(newItem, std::cout);
		}
		
		std::ostringstream out;
		printDAG(newItem, out);
		_resDAG = out.str();

		_ideotypeCost = newItemScore;
	}
}

void Solver::cross(const Genotype& genotype1, const Genotype& genotype2, 
	DpItemSet& newItemSet, bool updateFlag)
{
	static const int nLoci = _pData->getNumberOfLoci();
	static const DoubleMatrix& RM = _pData->getRM();
	static const double gamma = _pData->getGamma();
	static const unsigned long popMax = _pData->getPopMax();

	DpItem& item1 = _M->getItem(genotype1);
	DpItem& item2 = _M->getItem(genotype2);

	const GameteVector& gametes1 = item1.getGametes();
	const GameteVector& gametes2 = item2.getGametes();

	GenotypeSet ancestors;
	set_union(item1.getAncestors().begin(), item1.getAncestors().end(),
		item2.getAncestors().begin(), item2.getAncestors().end(),
		inserter(ancestors, ancestors.begin()));
	ancestors.insert(genotype1);
	ancestors.insert(genotype2);
	
	unsigned long cumPop = 0;
	unsigned long cumCross = 0;
	for (GenotypeSet::const_iterator it = ancestors.begin(); it != ancestors.end(); it++)
	{
		const DpItem& refItem = _M->getItem(*it);
		if (!refItem.isParent()) cumCross++;
		cumPop += refItem.getPop();
	}

	for (GameteVector::const_iterator it1 = gametes1.begin(); it1 != gametes1.end(); it1++)
	{
		for (GameteVector::const_iterator it2 = gametes2.begin(); it2 != gametes2.end(); it2++)
		{
			Genotype C(it1->_c, it2->_c);

			// kan er eigenlijk uit
			if (genotype1 == C || genotype2 == C || item1.hasAncestor(C) || item2.hasAncestor(C))
				continue;

			unsigned long pop = C.computePop(nLoci, RM, gamma, genotype1, genotype2);
			if (pop <= popMax)
			{
				DpItem newItem(C, &item1, &item2, pop, std::max(item1.getGen(), item2.getGen()) + 1, cumPop + pop, cumCross + 1, ancestors, updateFlag);
				//if (!_M->isPresent(C))
				//{
				//	DpItem dummyItem(C, NULL, NULL, INT_MAX, INT_MAX, INT_MAX, INT_MAX, GenotypeSet(), updateFlag);
				//	_M->updateItem(dummyItem);
				//}
				
				double newItemScore = getScore(newItem);
				if (newItemScore < _ideotypeCost)
				{
					if (!_M->isPresent(C) || newItemScore < getScore(_M->getItem(C)))
					//if (_currentGenotypes.find(C) == _currentGenotypes.end())
					{
						updateSet(newItem, newItemScore, newItemSet);
					}
					/*else 
					{
						assert(_M->isPresent(C));
						if (newItemScore < getScore(_M->getItem(C)))
						{
							updateSet(newItem, newItemScore, newItemSet);
						}
					}*/
				}
			}
		}
	}
}

double Solver::getScore(const DpItem& item) const
{
	return _pData->getCost(item.getCumPop(), item.getGen(), item.getCumCross());
}

void Solver::printDAG(const DpItem& item, std::ostream& out) const
{
	out << "digraph G {" << std::endl;
	item.printNodes(_M, _pData->getNumberOfLoci(), out);
	GenotypeSet genotypes;
	item.printEdges(out, _pData->getNumberOfLoci(), genotypes);
	out << "}" << std::endl;
}

std::string Solver::getMethodName() const
{
	return "DP";
}
