/*
 * genotypetable.cpp
 *
 *  Created on: 20-mar-2009
 *      Author: M. El-Kebir
 */

#include "genotypetable.h"
#include "data.h"

/* With m loci there are
 * 2^(2m-1) + 2^(m-1) different genotypes
 */

GenotypeTable* GenotypeTable::_pGenotypeTable = NULL;

GenotypeTable::GenotypeTable(unsigned int limit, int crossOverCount, bool restrictGametes)
	: _nLoci(Data::getInstance()->getNumberOfLoci())
	, _nGenotypes((1 << (2 * _nLoci - 1)) + (1 << (_nLoci - 1)))
	, _RM(Data::getInstance()->getRM())
	, _probLowerBound(Data::getInstance()->getProbLowerBound())
	, _ideotype(Data::getInstance()->getIdeotype())
	, _limit(limit)
	, _crossOverCount(crossOverCount)
	, _restrictGametes(restrictGametes)
	, _genotypeList()
{
	_pGenotypeTable = this;
}

GenotypeTable::~GenotypeTable()
{
	_pGenotypeTable = NULL;
}
