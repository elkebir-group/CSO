/*
 * dptable.cpp
 *
 *  Created on: 09-mar-2009
 *      Author: s030858
 */

#include "dptable.h"

/* With m loci there are
 * 2^(2m-1) + 2^(m-1) different genotypes
 */

DpTable::DpTable(int nLoci)
	: _nLoci(nLoci)
	, _nGenotypes((1 << (2 * _nLoci - 1)) + (1 << (_nLoci - 1)))
	, _useCount(0)
{
}

DpTable::~DpTable()
{
}
