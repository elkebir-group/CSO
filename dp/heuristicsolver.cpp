/*
 * heuristicsolver.cpp
 *
 *  Created on: 08-apr-2009
 *      Author: M. El-Kebir
 */

#include "heuristicsolver.h"
#include "heuristiccorrectloci.h"
#include "heuristiccorrectalleles.h"
#include "heuristicindividualcost.h"
#include "heuristiclinkage.h"
#include "heuristicpurelinkage.h"
#include "heuristiclargestsubgroupsize.h"
#include "heuristicsubgroupsizes.h"

HeuristicSolver::HeuristicSolver(const Data* pData, DpTable* pTable, int maxIterations, int verbosity, 
	HeuristicType heuristicType, int nGenotypesToSelect, int nGenotypesLinkageAnalysis,
	bool selfingRestriction)
	: Solver(pData, pTable, maxIterations, verbosity)
	, _heuristicType(heuristicType)
	, _nGenotypesToSelect(nGenotypesToSelect)
	, _linkageAnalysis(pData)
	, _nGenotypesLinkageAnalysis(nGenotypesLinkageAnalysis)
	, _selfingRestriction(selfingRestriction)
{
}

HeuristicSolver::~HeuristicSolver()
{
}

void HeuristicSolver::solve()
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

int HeuristicSolver::performStep(bool updateFlag)
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

void HeuristicSolver::applySortingHeuristic(DpItemVector& newItemVector, int target)
{
	switch (_heuristicType)
	{
	case HeuristicCorrectLociType:		// 1
		std::sort(newItemVector.begin(), newItemVector.end(), HeuristicCorrectLoci(target));
		break;
	case HeuristicCorrectAllelesType:	// 2
		std::sort(newItemVector.begin(), newItemVector.end(), HeuristicCorrectAlleles(target));
		break;
	case HeuristicLinkageType:	// 3
		std::sort(newItemVector.begin(), newItemVector.end(), HeuristicLinkage(target, HeuristicLinkageType));
		break;
	case HeuristicLinkageLociType:		// 4
		std::sort(newItemVector.begin(), newItemVector.end(), HeuristicLinkage(target, HeuristicCorrectLociType));
		break;
	case HeuristicLinkageAllelesType:	// 5
		std::sort(newItemVector.begin(), newItemVector.end(), HeuristicLinkage(target, HeuristicCorrectAllelesType));
		break;
	case HeuristicPureLinkageType:		// 10
		std::sort(newItemVector.begin(), newItemVector.end(), HeuristicPureLinkage(target, _linkageAnalysis));
		break;
	case HeuristicLargestSubGroupSizeType:	// 11
		std::sort(newItemVector.begin(), newItemVector.end(), HeuristicLargestSubGroupSize(target));
		break;
	case HeuristicSubGroupSizesType:	// 12
		std::sort(newItemVector.begin(), newItemVector.end(), HeuristicSubGroupSizes(target));
		break;
	}

	if (_verbosity > 1)
	{
		int counter = 1;
		for (DpItemVector::const_iterator it = newItemVector.begin(); it != newItemVector.end(); it++, counter++)
		{
			std::cout << "// " << counter << "\t";
			it->printItem(_pData->getNumberOfLoci(), std::cout, false);
			switch (_heuristicType)
			{
			case HeuristicCorrectLociType:		// 1
				std::cout << "\t" << HeuristicCorrectLoci(target).getCost(*it) << std::endl;
				break;
			case HeuristicCorrectAllelesType:	// 2
				std::cout << "\t" << HeuristicCorrectAlleles(target).getCost(*it) << std::endl;
				break;
			case HeuristicLinkageType:	// 3
				std::cout << "\t" << HeuristicLinkage(target, HeuristicLinkageType).getCost(*it) << std::endl;
				break;
			case HeuristicLinkageLociType:		// 4
				std::cout << "\t" << HeuristicLinkage(target, HeuristicCorrectLociType).getCost(*it) << std::endl;
				break;
			case HeuristicLinkageAllelesType:	// 5
				std::cout << "\t" << HeuristicLinkage(target, HeuristicCorrectAllelesType).getCost(*it) << std::endl;
				break;
			case HeuristicPureLinkageType:		// 10
				std::cout << "\t" << HeuristicPureLinkage(target, _linkageAnalysis).getCost(*it) << std::endl;
				break;
			case HeuristicLargestSubGroupSizeType:		// 11
				std::cout << "\t" << HeuristicLargestSubGroupSize(target).getCost(*it) << std::endl;
				break;
			case HeuristicSubGroupSizesType:		// 12
				std::cout << "\t" << HeuristicSubGroupSizes(target).getCost(*it) << std::endl;
				break;
			}
		}
	}
}

void HeuristicSolver::applyLinkageAnalysisUnlinked(const DpItemVector& newItemVector, const LociPairList& lociList, 
	const int targetChromosome, DpItemSet& itemsToKeep) const
{
	for (LociPairList::const_iterator lociIt = lociList.begin(); lociIt != lociList.end(); lociIt++)
	{
		/* Search for a genotype in which the loci are:
		 * - Linked
		 * - WeaklyLinked
		 */

		int count = 0;
		int foundLinked = 0;
		int foundWeaklyLinked = 0;
		for (DpItemVector::const_iterator it = newItemVector.begin();
			it != newItemVector.end() && !(foundLinked == _nGenotypesLinkageAnalysis && foundWeaklyLinked == _nGenotypesLinkageAnalysis); it++)
		{
			switch (it->getLinkage(lociIt->first, lociIt->second, targetChromosome))
			{
			case Linked:
				if (foundLinked < _nGenotypesLinkageAnalysis)
				{
					foundLinked++;
					if (count >= _nGenotypesToSelect)
						itemsToKeep.insert(*it);
				}
				break;
			case WeaklyLinked:
				if (foundWeaklyLinked < _nGenotypesLinkageAnalysis)
				{
					foundWeaklyLinked++;
					if (count >= _nGenotypesToSelect)
						itemsToKeep.insert(*it);
				}
				break;
			case Unlinked:
				// do nothing
				break;
			}

			count++;
		}
	}
}

void HeuristicSolver::applyLinkageAnalysisWeaklyLinked(const DpItemVector& newItemVector, const LociPairList& lociList, 
	const int targetChromosome, DpItemSet& itemsToKeep) const
{
	for (LociPairList::const_iterator lociIt = lociList.begin(); lociIt != lociList.end(); lociIt++)
	{
		/* Search for a genotype in which the loci are:
		 * - Linked
		 */

		int count = 0;
		int foundLinked = 0;
		for (DpItemVector::const_iterator it = newItemVector.begin();
			it != newItemVector.end() && foundLinked != _nGenotypesLinkageAnalysis; it++)
		{
			switch (it->getLinkage(lociIt->first, lociIt->second, targetChromosome))
			{
			case Linked:
				if (foundLinked < _nGenotypesLinkageAnalysis)
				{
					foundLinked++;
					if (count >= _nGenotypesToSelect)
						itemsToKeep.insert(*it);
				}
				break;
			case WeaklyLinked:
			case Unlinked:
				// do nothing
				break;
			}

			count++;
		}
	}
}

std::string HeuristicSolver::getMethodName() const
{
	char buf[1024];

	if (_nGenotypesLinkageAnalysis)
	{
		if (_M->getRestrictGametes())
		{
			sprintf(buf, "DP%d-%d-%dg", (int) _heuristicType, _nGenotypesToSelect, _nGenotypesLinkageAnalysis);
		}
		else
		{
			sprintf(buf, "DP%d-%d-%d", (int) _heuristicType, _nGenotypesToSelect, _nGenotypesLinkageAnalysis);
		}
	}
	else
	{
		if (_M->getRestrictGametes())
		{
			sprintf(buf, "DP%d-%dg", (int) _heuristicType, _nGenotypesToSelect);
		}
		else
		{
			sprintf(buf, "DP%d-%d", (int) _heuristicType, _nGenotypesToSelect);
		}
	}

	return std::string(buf);
}

int HeuristicSolver::applyInitialParentSelfing()
{
	DpItemSet newItemSet;

	const GenotypeSet& parents = Data::getInstance()->getParents();
	const unsigned long popMax = Data::getInstance()->getPopMax();
	const DoubleMatrix& RM = Data::getInstance()->getRM();
	const double gamma = Data::getInstance()->getGamma();
	const int nLoci = Data::getInstance()->getNumberOfLoci();

	if (!Data::getInstance()->getAllowSelfing())
		return 0;	// selfing is not allowed

	for (GenotypeSet::const_iterator it = parents.begin(); it != parents.end(); it++)
	{
		if (it->isHomozygous())
			continue;	// no need for selfing, parent is already homozygous

		DpItem& item = _M->getItem(*it);

		Genotype C(it->getC0(), it->getC0());
		Genotype D(it->getC1(), it->getC1());

		if (parents.find(C) == parents.end())
		{
			unsigned long pop = C.computePop(nLoci, RM, gamma, *it, *it);
			if (0 <= pop && pop <= popMax)
			{
				GenotypeSet ancestors;
				ancestors.insert(*it);

				DpItem newItem(C, &item, &item, pop, 
					item.getGen() + 1, 
					item.getCumPop() + pop, 
					item.getCumCross() + 1, ancestors, false);
				
				updateSet(newItem, getScore(newItem), newItemSet);
			}
		}

		if (parents.find(D) == parents.end())
		{
			unsigned long pop = D.computePop(nLoci, RM, gamma, *it, *it);
			if (0 <= pop && pop <= popMax)
			{
				GenotypeSet ancestors;
				ancestors.insert(*it);

				DpItem newItem(D, &item, &item, pop, 
					item.getGen() + 1, 
					item.getCumPop() + pop, 
					item.getCumCross() + 1, ancestors, false);
				
				updateSet(newItem, getScore(newItem), newItemSet);
			}
		}
	}

	_currentGenotypes.insert(newItemSet.begin(), newItemSet.end());

	for (DpItemSet::iterator it = newItemSet.begin(); it != newItemSet.end(); it++)
	{
		_M->updateItem(*it).updateAttributesFast(_M, true);
	}

	if (_verbosity > 0)
	{
		std::cout << "// Initial parent selfing yielded " << newItemSet.size() << " genotypes." << std::endl;
	}

	return (int) newItemSet.size();
}

void HeuristicSolver::cross(const Genotype& genotype1, const Genotype& genotype2, 
	DpItemSet& newItemSet, bool updateFlag)
{
	static const int nLoci = _pData->getNumberOfLoci();
	static const DoubleMatrix& RM = _pData->getRM();
	static const double gamma = _pData->getGamma();
	static const unsigned long popMax = _pData->getPopMax();
	static const bool allowSelfing = _pData->getAllowSelfing();

	DpItem& item1 = _M->getItem(genotype1);
	DpItem& item2 = _M->getItem(genotype2);

	const GameteVector& gametes1 = item1.getGametes();
	const GameteVector& gametes2 = item2.getGametes();

	/* generate ancestor set */
	GenotypeSet ancestors;
	set_union(item1.getAncestors().begin(), item1.getAncestors().end(),
		item2.getAncestors().begin(), item2.getAncestors().end(),
		inserter(ancestors, ancestors.begin()));
	ancestors.insert(genotype1);
	ancestors.insert(genotype2);
	
	/* set cumPop and cumCross */
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

			// no intermediate selfing
			if (_selfingRestriction && item1 == item2 && C != _pData->getIdeotype())
				continue;

			// if C is homozygous and not obtained via a selfing then it is not optimal
			if (allowSelfing && C.isHomozygous() && item1 != item2)
				continue;

			unsigned long pop = C.computePop(nLoci, RM, gamma, genotype1, genotype2);
			if (pop <= popMax)
			{
				DpItem newItem(C, &item1, &item2, 
					pop, 
					std::max(item1.getGen(), item2.getGen()) + 1, 
					cumPop + pop, 
					cumCross + 1, 
					ancestors, updateFlag);
				
				double newItemScore = getScore(newItem);
				if (newItemScore < _ideotypeCost)
				{
					if (!_M->isPresent(C) || newItemScore < getScore(_M->getItem(C)))
					{
						updateSet(newItem, newItemScore, newItemSet);
					}
				}
			}
		}
	}
}

DpItemSet HeuristicSolver::applyHeuristic(const DpItemSet& newItemSet, int target)
{
	DpItemVector newItemVector(newItemSet.begin(), newItemSet.end());
	applySortingHeuristic(newItemVector, target);

	DpItemSet itemsToKeep;
	if (_nGenotypesLinkageAnalysis)
	{
		const LociPairList& unlinkedLoci0 = _linkageAnalysis.getUnlinkedLoci(target);
		applyLinkageAnalysisUnlinked(newItemVector, unlinkedLoci0, target, itemsToKeep);

		if (_verbosity > 0)
		{
			std::cout << "// Number of selected genotypes due to linkage analysis: " << itemsToKeep.size() << std::endl;
		}
	}

	if ((int) newItemVector.size() > _nGenotypesToSelect)
		newItemVector.erase(newItemVector.begin() + _nGenotypesToSelect, newItemVector.end());

	itemsToKeep.insert(newItemVector.begin(), newItemVector.end());

	/* make sure ideotype is included in itemsToKeep, if it is present in newItemSet */
	const Genotype& ideotype = _pData->getIdeotype();
	DpItem ideotypeDummyItem = DpItem(ideotype, false);

	DpItemSet::const_iterator ideotypeIt = newItemSet.find(ideotypeDummyItem);
	if (ideotypeIt != newItemSet.end() && itemsToKeep.find(ideotypeDummyItem) == itemsToKeep.end())
	{
		if (_verbosity > 0)
		{
			std::cout << "// Manually added the ideotype as it was not selected by the heuristic." << std::endl;
		}
		itemsToKeep.insert(*ideotypeIt);
	}

	return itemsToKeep;
}