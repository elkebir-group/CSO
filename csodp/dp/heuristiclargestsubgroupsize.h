/*
 * heuristiclargestsubgroupsize.h
 *
 *  Created on: 09-may-2009
 *      Author: M. El-Kebir
 */

#ifndef HEURISTICLARGESTSUBGROUPSIZE_H
#define HEURISTICLARGESTSUBGROUPSIZE_H

#include <functional>
#include <algorithm>
#include "cso.h"
#include "data.h"
#include "genotype.h"
#include "dpitem.h"
#include "heuristicbase.h"

typedef std::pair<int, int> GameteGroupSizePair;

class HeuristicLargestSubGroupSize : public HeuristicBase
{
private:
	int getNumberOfCorrectAlleles(const DpItem& item) const;
	GameteGroupSizePair getInternalCost(const DpItem& item) const;

public:
	HeuristicLargestSubGroupSize(int target);
	bool operator ()(const DpItem& item1, const DpItem& item2) const;
	int getCost(const DpItem& item) const;
};

inline HeuristicLargestSubGroupSize::HeuristicLargestSubGroupSize(int target)
	: HeuristicBase(target)
{
}

inline bool HeuristicLargestSubGroupSize::operator ()(const DpItem& item1, const DpItem& item2) const
{
	GameteGroupSizePair cost1 = getInternalCost(item1);
	GameteGroupSizePair cost2 = getInternalCost(item2);
	
	if (cost1.second < cost2.second)
		return false;
	else if (cost1.second > cost2.second)
		return true;
	else
	{
		//int costAlleles1 = getNumberOfCorrectAlleles(item1);
		//int costAlleles2 = getNumberOfCorrectAlleles(item2);

		int costAlleles1 = _nLoci - numberOfDifferences(_nLoci, cost1.first, _target);
		int costAlleles2 = _nLoci - numberOfDifferences(_nLoci, cost2.first, _target); 

		if (costAlleles1 < costAlleles2)
			return false;
		else if (costAlleles1 > costAlleles2)
			return true;
		else
		{
			// cumPop + pop werkt beter ;)
			//return item1.computeProb(nLoci, RM, cost1.first) > item2.computeProb(nLoci, RM, cost2.first);
			return (item1.getCumPop() + probToPop(item1.computeProb(_nLoci, _RM, cost1.first), _gamma)) < 
				(item2.getCumPop() + probToPop(item2.computeProb(_nLoci, _RM, cost2.first), _gamma));
		}
	}
}

inline GameteGroupSizePair HeuristicLargestSubGroupSize::getInternalCost(const DpItem& item) const
{
	int gamete = item.getMaskAtMostOneCrossOverGroup(_nLoci, _target);
	return GameteGroupSizePair(gamete, largestSubGroupSize(_nLoci, gamete, _target));
}

inline int HeuristicLargestSubGroupSize::getCost(const DpItem& item) const
{
	return getInternalCost(item).second;
}

inline int HeuristicLargestSubGroupSize::getNumberOfCorrectAlleles(const DpItem& item) const
{
	return _nLoci - std::min(numberOfDifferences(_nLoci, item.getC0(), _target), 
		numberOfDifferences(_nLoci, item.getC1(), _target));
}

#endif /* HEURISTICLARGESTSUBGROUPSIZE_H */