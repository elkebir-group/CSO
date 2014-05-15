/*
 * dptable.cpp
 *
 *  Created on: 09-mar-2009
 *      Author: M. El-Kebir
 */

#include "dptable.h"
#include "data.h"

/* With m loci there are
 * 2^(2m-1) + 2^(m-1) different genotypes
 */

DpTable::DpTable(int crossOverCount, bool restrictGametes)
  : _nLoci(Data::getInstance()->getNumberOfLoci())
  , _nGenotypes((1 << (2 * _nLoci - 1)) + (1 << (_nLoci - 1)))
  , _RM(Data::getInstance()->getRM())
  , _probLowerBound(Data::getInstance()->getProbLowerBound())
  , _ideotype(Data::getInstance()->getIdeotype())
  , _crossOverCount(crossOverCount)
  , _restrictGametes(restrictGametes)
  , _useCount(0)
{
}

DpTable::~DpTable()
{
}
