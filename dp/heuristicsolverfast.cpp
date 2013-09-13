/*
 * heuristicsolverfast.cpp
 *
 *  Created on: 09-apr-2009
 *      Author: M. El-Kebir
 */

#include "heuristicsolverfast.h"

HeuristicSolverFast::HeuristicSolverFast(const Data* pData, DpTable* pTable, 
		int maxIterations, int verbosity, 
		HeuristicType heuristicType, int nGenotypesToSelect, 
		int nGenotypesLinkageAnalysis, bool selfingRestriction)
	: HeuristicSolver(pData, pTable, 
		maxIterations, verbosity, 
		heuristicType, nGenotypesToSelect, 
		nGenotypesLinkageAnalysis, selfingRestriction)
{
}

HeuristicSolverFast::~HeuristicSolverFast()
{
}

void HeuristicSolverFast::solve()
{
	applyInitialParentSelfing();

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

int HeuristicSolverFast::performStep(bool updateFlag)
{
	const Genotype& ideotype = _pData->getIdeotype();
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
				/* heuristic: do not self a non-parent or a non-final genotype */
				if (ideotype.isHomozygous() && it2 == it1)
				{
					const DpItem& item = _M->getItem(*it1);
					if (!item.isParent() &&
						numberOfDifferences(_pData->getNumberOfLoci(), item.getC0() | item.getC1(), 
							_pData->getIdeotype().getC0()) != 0)
					{
						continue;
					}
				}
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

	/* perform crossings among current and old individuals, these two sets are disjoint */
	for (GenotypeSet::iterator it1 = _currentGenotypes.begin(); it1 != _currentGenotypes.end(); it1++)
	{
		for (GenotypeSet::iterator it2 = _oldGenotypes.begin(); it2 != _oldGenotypes.end(); it2++)
		{
			cross(*it1, *it2, newItemSet, updateFlag); // O(4^m)
		}
	}

	if ((erasedCount = cleanUp(newItemSet)))
	{
		if (_verbosity > 0)
			std::cout << "// Erased " << erasedCount << " new sub-optimal genotypes." << std::endl;
	}

	DpItemSet itemsToKeep = applyHeuristic(newItemSet, ideotype.getC0());;
	if (!ideotype.isHomozygous())
	{
		DpItemSet itemsToKeep2 = applyHeuristic(newItemSet, ideotype.getC1());
		itemsToKeep.insert(itemsToKeep2.begin(), itemsToKeep2.end());
	}

	int cheaper = 0;
	for (DpItemSet::iterator it = itemsToKeep.begin(); it != itemsToKeep.end(); it++)
	{
		cheaper += _currentGenotypes.erase(*it);
		cheaper += _oldGenotypes.erase(*it);
		
		_M->updateItem(*it);
	}

	if (_verbosity > 0)
	{
		std::cout << "// Number of new individuals: " << newItemSet.size() - cheaper << std::endl;
		std::cout << "// Number of selected cheaper individuals: " << cheaper << std::endl;
	}

	_oldGenotypes.insert(_currentGenotypes.begin(), _currentGenotypes.end());
	if ((erasedCount = cleanUp(_oldGenotypes)))
	{
		if (_verbosity > 0)
			std::cout << "// Erased " << erasedCount << " old sub-optimal genotypes." << std::endl;
	}

	_currentGenotypes.clear();
	_currentGenotypes.insert(itemsToKeep.begin(), itemsToKeep.end());

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
