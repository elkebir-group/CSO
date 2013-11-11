/*
 * heuristiclinkage.h
 *
 *  Created on: 07-apr-2009
 *      Author: M. El-Kebir
 */

#ifndef HEURISTICLINKAGE_H
#define HEURISTICLINKAGE_H

#include <functional>
#include <algorithm>
#include "cso.h"
#include "data.h"
#include "genotype.h"
#include "dpitem.h"
#include "heuristicbase.h"
#include "heuristiccorrectalleles.h"
#include "heuristiccorrectloci.h"

class HeuristicLinkage : public HeuristicBase
{
private:
	HeuristicType _heuristicType;
	const double _penaltyCost;

public:
	HeuristicLinkage(int target, HeuristicType heuristicType);
	bool operator ()(const DpItem& item1, const DpItem& item2) const;
	double getCost(const DpItem& item) const;
};

inline HeuristicLinkage::HeuristicLinkage(int target, HeuristicType heuristicType)
	: HeuristicBase(target)
	, _heuristicType(heuristicType)
	, _penaltyCost(std::max(Data::getInstance()->getCostCross(), Data::getInstance()->getCostGen()))
{
}

inline double HeuristicLinkage::getCost(const DpItem& item) const
{
	int c0 = item.getC0();
	int c1 = item.getC1();

	int penalty = 0;
	for (int i = 0; i < _nLoci; i++)
	{
        if (GET_BIT(_nLoci, c0, i) != GET_BIT(_nLoci, _target, i) &&
            GET_BIT(_nLoci, c1, i) != GET_BIT(_nLoci, _target, i))
		{
			penalty++;
		}
	}

	//int mask = item.getMaskAtMostOneCrossOverGroup(nLoci, target);
	int mask = item.getMask(_nLoci, _target);
	
	double nominator = probToPop(item.computeProb(_nLoci, _RM, mask), _gamma) + 
		penalty * _penaltyCost + item.getPop() / 10.0;

	switch (_heuristicType)
	{
	case HeuristicCorrectLociType:
		return nominator / (HeuristicCorrectLoci(_target).getCost(item) + DBL_EPSILON);
	case HeuristicCorrectAllelesType:
		return nominator / (HeuristicCorrectAlleles(_target).getCost(item) + DBL_EPSILON);
	case HeuristicLinkageType:
		return nominator;
	default:
		assert(false);
		return DBL_MAX;
	}

}

inline bool HeuristicLinkage::operator ()(const DpItem& item1, const DpItem& item2) const
{
	double cost1 = getCost(item1);
	double cost2 = getCost(item2);
	
	if (cost1 > cost2)
		return false;
	else if (cost1 < cost2)
		return true;
	else
		return item1.compare(item2) < 0;
}

#endif /* HEURISTICLINKAGE_H */
