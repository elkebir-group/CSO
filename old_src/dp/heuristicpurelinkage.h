/*
 * heuristicpurelinkage.h
 *
 *  Created on: 05-may-2009
 *      Author: M. El-Kebir
 */

#ifndef HEURISTICPURELINKAGE_H
#define HEURISTICPURELINKAGE_H

#include <functional>
#include <algorithm>
#include "cso.h"
#include "data.h"
#include "genotype.h"
#include "dpitem.h"
#include "linkageanalysis.h"
#include "heuristicbase.h"
#include "heuristiclargestsubgroupsize.h"

class HeuristicPureLinkage : public HeuristicBase
{
private:
	const LinkageAnalysis& _linkageAnalysis;
	static const int _costMatrix[3][3];

public:
	HeuristicPureLinkage(int target, const LinkageAnalysis& linkageAnalysis);
	bool operator ()(const DpItem& item1, const DpItem& item2) const;
	double getCost(const DpItem& item) const;
};

inline bool HeuristicPureLinkage::operator ()(const DpItem& item1, const DpItem& item2) const
{
	double cost1 = getCost(item1);
	double cost2 = getCost(item2);
	
	if (cost1 > cost2)
		return false;
	else if (cost1 < cost2)
		return true;
	else
		return HeuristicLargestSubGroupSize(_target)(item1, item2);
}

#endif
