/*
 * heuristicindividualcost.h
 *
 *  Created on: 03-apr-2009
 *      Author: M. El-Kebir
 */

#ifndef HEURISTICINDIVIDUALCOST_H
#define HEURISTICINDIVIDUALCOST_H

#include <functional>
#include "cso.h"
#include "data.h"
#include "dpitem.h"

class HeuristicIndividualCost : public std::binary_function<DpItem, DpItem, bool>
{
public:
  bool operator ()(const DpItem& item1, const DpItem& item2) const;
  double getCost(const DpItem& item) const;
};

inline bool HeuristicIndividualCost::operator ()(const DpItem& item1, const DpItem& item2) const
{
  double cost1 = getCost(item1);
  double cost2 = getCost(item2);

  if (cost1 < cost2)
    return true;
  else if (cost1 > cost2)
    return false;
  else
    return item1.compare(item2) < 0;
}

inline double HeuristicIndividualCost::getCost(const DpItem& item) const
{
  return Data::getInstance()->getCost(item.getCumPop(), item.getGen(), item.getCumCross());
}

#endif /* HEURISTICINDIVIDUALCOST_H */
