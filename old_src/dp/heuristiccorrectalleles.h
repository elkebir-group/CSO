/*
 * heuristiccorrectalleles.h
 *
 *  Created on: 07-apr-2009
 *      Author: M. El-Kebir
 */

#ifndef HEURISTICCORRECTALLELES_H
#define HEURISTICCORRECTALLELES_H

#include <functional>
#include <algorithm>
#include "cso.h"
#include "data.h"
#include "genotype.h"
#include "dpitem.h"
#include "heuristicbase.h"
#include "heuristiclargestsubgroupsize.h"

class HeuristicCorrectAlleles : public HeuristicBase
{
public:
  HeuristicCorrectAlleles(int target);
  bool operator ()(const DpItem& item1, const DpItem& item2) const;
  int getCost(const DpItem& item) const;
};

inline HeuristicCorrectAlleles::HeuristicCorrectAlleles(int target)
  : HeuristicBase(target)
{
}

inline bool HeuristicCorrectAlleles::operator ()(const DpItem& item1, const DpItem& item2) const
{
  double cost1 = getCost(item1);
  double cost2 = getCost(item2);

  if (cost1 < cost2)
    return false;
  else if (cost1 > cost2)
    return true;
  else
    return HeuristicLargestSubGroupSize(_target)(item1, item2);
}

inline int HeuristicCorrectAlleles::getCost(const DpItem& item) const
{
  return _nLoci - std::min(numberOfDifferences(_nLoci, item.getC0(), _target),
    numberOfDifferences(_nLoci, item.getC1(), _target));
}

#endif
