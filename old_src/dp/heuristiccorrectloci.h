/*
 * heuristiccorrectloci.h
 *
 *  Created on: 07-apr-2009
 *      Author: M. El-Kebir
 */

#ifndef HEURISTICCORRECTLOCI_H
#define HEURISTICCORRECTLOCI_H

#include <functional>
#include "cso.h"
#include "data.h"
#include "genotype.h"
#include "dpitem.h"
#include "heuristicbase.h"
#include "heuristiccorrectalleles.h"
#include "heuristiclargestsubgroupsize.h"

class HeuristicCorrectLoci : public HeuristicBase
{
public:
	HeuristicCorrectLoci(int target);
	bool operator ()(const DpItem& item1, const DpItem& item2) const;
	int getCost(const DpItem& item) const;
};

inline HeuristicCorrectLoci::HeuristicCorrectLoci(int target)
	: HeuristicBase(target)
{
}

inline bool HeuristicCorrectLoci::operator ()(const DpItem& item1, const DpItem& item2) const
{
	double cost1 = getCost(item1);
	double cost2 = getCost(item2);
	
	if (cost1 < cost2)
		return false;
	else if (cost1 > cost2)
		return true;
	else
		return HeuristicCorrectAlleles(_target)(item1, item2);
	/*{
		double costAlleles1 = HeuristicCorrectAlleles::getCost(item1);
		double costAlleles2 = HeuristicCorrectAlleles::getCost(item2);

		if (costAlleles1 < costAlleles2)
			return false;
		else if (costAlleles1 > costAlleles2)
			return true;
		else
		{
			double costLargestSubGroup1 = HeuristicLargestSubGroupSize::getCost(item1);
			double costLargestSubGroup2 = HeuristicLargestSubGroupSize::getCost(item2);

			return costLargestSubGroup1 < costLargestSubGroup2;
		}
	}*/
}

inline int HeuristicCorrectLoci::getCost(const DpItem& item) const
{
	int nCorrectLoci = 0;
	for (int i = 0; i < _nLoci; i++)
	{
		if (item(0, i) == GET_BIT(_target, i) || item(1, i) == GET_BIT(_target, i))
			nCorrectLoci++;
	}

	return nCorrectLoci;
}

#endif /* HEURISTICCORRECTLOCI_H */
