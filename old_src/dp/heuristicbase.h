/*
 * heuristicbase.h
 *
 *  Created on: 13-may-2009
 *      Author: M. El-Kebir
 */

#ifndef HEURISTICBASE_H
#define HEURISTICBASE_H

#include <functional>
#include <algorithm>
#include "cso.h"
#include "data.h"
#include "genotype.h"
#include "dpitem.h"

class HeuristicBase : public std::binary_function<DpItem, DpItem, bool>
{
protected:
  const int _target;
  const int _nLoci;
  const DoubleMatrix& _RM;
  const double _gamma;

public:
  HeuristicBase(int target);
  virtual bool operator ()(const DpItem& item1, const DpItem& item2) const = 0;
};

inline HeuristicBase::HeuristicBase(int target)
  : _target(target)
  , _nLoci(Data::getInstance()->getNumberOfLoci())
  , _RM(Data::getInstance()->getRM())
  , _gamma(Data::getInstance()->getGamma())
{
}

#endif /* HEURISTICBASE_H */
