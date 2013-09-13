/*
 * heuristicsubgroupsize.h
 *
 *  Created on: 09-may-2009
 *      Author: M. El-Kebir
 */

#ifndef HEURISTICSUBGROUPSIZE_H
#define HEURISTICSUBGROUPSIZE_H

#include <functional>
#include <algorithm>
#include "cso.h"
#include "data.h"
#include "genotype.h"
#include "dpitem.h"
#include "heuristiccorrectalleles.h"

typedef std::pair<int, int> GameteGroupSizePair;

class HeuristicSubGroupSize : public std::binary_function<DpItem, DpItem, bool>
{
private:
	static GameteGroupSizePair compute(int nLoci, const DpItem& item, int target);

public:
	bool operator ()(const DpItem& item1, const DpItem& item2) const;
	static double getCost(const DpItem& item);
	static GameteGroupSizePair getInternalCost(const DpItem& item);
};

inline bool HeuristicSubGroupSize::operator ()(const DpItem& item1, const DpItem& item2) const
{
	static const int nLoci = Data::getInstance()->getNumberOfLoci();
	static const DoubleMatrix& RM = Data::getInstance()->getRM();
	static const double gamma = Data::getInstance()->getGamma();

	GameteGroupSizePair cost1 = getInternalCost(item1);
	GameteGroupSizePair cost2 = getInternalCost(item2);
	
	if (cost1.second < cost2.second)
		return false;
	else if (cost1.second > cost2.second)
		return true;
	else
	{
		double costAlleles1 = HeuristicCorrectAlleles::getCost(item1);
		double costAlleles2 = HeuristicCorrectAlleles::getCost(item2);

		if (costAlleles1 < costAlleles2)
			return false;
		else if (costAlleles1 > costAlleles2)
			return true;
		else
		{
			// cumPop + pop werkt beter ;)
			//return item1.computeProb(nLoci, RM, cost1.first) > item2.computeProb(nLoci, RM, cost2.first);
			return (item1.getCumPop() + probToPop(item1.computeProb(nLoci, RM, cost1.first), gamma)) < 
				(item2.getCumPop() + probToPop(item2.computeProb(nLoci, RM, cost2.first), gamma));
		}
	}
}

inline GameteGroupSizePair HeuristicSubGroupSize::compute(int nLoci, const DpItem& item, int target)
{
	int gamete = item.getMaskAtMostOneCrossOverGroup(nLoci, target);
	return GameteGroupSizePair(gamete, largestSubGroupSize(nLoci, gamete, target));
}

inline GameteGroupSizePair HeuristicSubGroupSize::getInternalCost(const DpItem& item)
{
	static const Genotype& ideotype = Data::getInstance()->getIdeotype();
	static const int nLoci = Data::getInstance()->getNumberOfLoci();

	GameteGroupSizePair pair0 = compute(nLoci, item, ideotype.getC0());
	GameteGroupSizePair pair1 = compute(nLoci, item, ideotype.getC1());
	
	if (pair0 > pair1)
	{
		return pair0;
	}
	else
	{
		return pair1;
	}
}

inline double HeuristicSubGroupSize::getCost(const DpItem& item)
{
	return getInternalCost(item).second;
}


#endif /* HEURISTICSUBGROUPSIZE_H */