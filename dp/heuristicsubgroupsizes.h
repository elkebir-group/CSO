/*
 * heuristicsubgroupsizes.h
 *
 * TODO: aanpassen zodat deze werkt met heterozygote clusters...
 *
 *  Created on: 09-may-2009
 *      Author: M. El-Kebir
 */

#ifndef HEURISTICSUBGROUPSIZES_H
#define HEURISTICSUBGROUPSIZES_H

#include <functional>
#include <algorithm>
#include "cso.h"
#include "data.h"
#include "genotype.h"
#include "dpitem.h"
#include "heuristicbase.h"
#include "heuristiccorrectalleles.h"

typedef std::pair<int, int> GameteScorePair;

class HeuristicSubGroupSizes : public HeuristicBase
{
private:
	GameteScorePair compute(int gamete) const;
	GameteScorePair getInternalCost(const DpItem& item) const;

public:
	HeuristicSubGroupSizes(int target);
	bool operator ()(const DpItem& item1, const DpItem& item2) const;
	int getCost(const DpItem& item) const;	
};

inline HeuristicSubGroupSizes::HeuristicSubGroupSizes(int target)
	: HeuristicBase(target)
{
}

inline bool HeuristicSubGroupSizes::operator ()(const DpItem& item1, const DpItem& item2) const
{
	/*GameteScorePair gameteScorePair1 = getInternalCost(item1);
	GameteScorePair gameteScorePair2 = getInternalCost(item2);

	if (gameteScorePair1.second < gameteScorePair2.second)
		return false;
	else if (gameteScorePair1.second > gameteScorePair2.second)
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
			return (item1.getCumPop() + probToPop(item1.computeProb(nLoci, RM, gameteScorePair1.first), gamma)) < 
				(item2.getCumPop() + probToPop(item2.computeProb(nLoci, RM, gameteScorePair2.first), gamma));
		}
	}*/
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

inline GameteScorePair HeuristicSubGroupSizes::compute(int gamete) const
{
	int score = 0;
	int maxGroupSize = 0;
	for (int i = 0; i < _nLoci; i++)
	{
		if (GET_BIT(gamete, i) == GET_BIT(_target, i))
		{
			maxGroupSize++;
		}
		else if (maxGroupSize != 0)
		{
			score += (int) pow((double) 10, maxGroupSize - 1);
			maxGroupSize = 0;
		}
	}
	
	if (maxGroupSize != 0)
	{
		score += (int) pow((double) 10, maxGroupSize - 1);
		maxGroupSize = 0;
	}
	return GameteScorePair(gamete, score);
}

inline GameteScorePair HeuristicSubGroupSizes::getInternalCost(const DpItem& item) const
{
	int c0 = item.getC0();
	int c1 = item.getC1();

	std::vector<int> homozygousLoci;
	std::vector<int> heterozygousLoci;

	std::vector<int> crossOverPoints;
	chromosomeCompare(_nLoci, c0, c1, homozygousLoci, heterozygousLoci);
	for (int i = 0; i < (int) heterozygousLoci.size(); i++)
	{
		int val = heterozygousLoci[i];
		if (val == 0)
		{
			crossOverPoints.push_back(val);
		}
		else if (val == _nLoci - 1)
		{
			if (crossOverPoints.size() == 0 || crossOverPoints.back() != val - 1)
				crossOverPoints.push_back(val - 1);
		}
		else
		{
			if (crossOverPoints.size() == 0 || crossOverPoints.back() != val - 1)
				crossOverPoints.push_back(val - 1);
			crossOverPoints.push_back(val);
		}
	}

	GameteScorePair gameteScorePair0 = compute(c0);
	GameteScorePair gameteScorePair1 = compute(c1);
	
	GameteScorePair maxGameteScorePair;
	if (gameteScorePair0.second > gameteScorePair1.second)
	{
		maxGameteScorePair = gameteScorePair0;
	}
	else
	{
		maxGameteScorePair = gameteScorePair1;
	}
	
	for (int i = 0; i < (int) crossOverPoints.size(); i++)
	{
		int mask = GENERATE_AND_MASK(0, crossOverPoints[i] + 1);
		int compMask = (~mask) & (GENERATE_AND_MASK(0, _nLoci));

		int gamete0 = (c0 & mask) | (c1 & compMask);
		int gamete1 = (c1 & mask) | (c0 & compMask);

		gameteScorePair0 = compute(gamete0);
		gameteScorePair1 = compute(gamete1);

		if (gameteScorePair0.second > maxGameteScorePair.second)
		{
			maxGameteScorePair = gameteScorePair0;
		}
		else if (gameteScorePair1.second > maxGameteScorePair.second)
		{
			maxGameteScorePair = gameteScorePair1;
		}
	}

	return maxGameteScorePair;
}

inline int HeuristicSubGroupSizes::getCost(const DpItem& item) const
{
	return getInternalCost(item).second;
}

#endif /* HEURISTICSUBGROUPSIZES_H */
