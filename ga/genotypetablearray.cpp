/*
 * genotypetablearray.cpp
 *
 *  Created on: 20-mar-2009
 *      Author: M. El-Kebir
 */

#include "genotypetablearray.h"
#include <string.h>

GenotypeTableArray::GenotypeTableArray(unsigned int limit, int crossOverCount, bool restrictGametes)
	: GenotypeTable(limit, crossOverCount, restrictGametes)
	, _A(new GenotypeGamete*[_nGenotypes])
{
	memset(_A, NULL, sizeof(GenotypeGamete*) * _nGenotypes);
}

GenotypeTableArray::~GenotypeTableArray()
{
	for (int i = 0; i < _nGenotypes; i++)
	{
		delete _A[i];
	}

	delete[] _A;
}

const GenotypeGamete& GenotypeTableArray::getGenotype(const Genotype& genotype) const
{
	GenotypeGamete* pResult = _A[getIndex(genotype)];
	if (!pResult)
	{
		pResult = new GenotypeGamete(genotype.getC0(), genotype.getC1());

		if (_crossOverCount > 0)
			pResult->computeGametes(_nLoci, _RM, _ideotype, 
				_probLowerBound, _restrictGametes, _crossOverCount);
		else
			pResult->computeGametes(_nLoci, _RM, _probLowerBound);

		pResult->computeGametesCumulative();
		_A[getIndex(genotype)] = pResult;
		
		_genotypeList.push_back(genotype);
	}

	if (_genotypeList.size() == _limit + 1)
	{
		Genotype oldestGenotype = _genotypeList.front();
		_genotypeList.pop_front();
		delete _A[getIndex(oldestGenotype)];
		_A[getIndex(oldestGenotype)] = NULL;
	}

	assert(_genotypeList.size() <= _limit);

	return *pResult;
}

const GenotypeGamete& GenotypeTableArray::getGenotype(int c0, int c1) const
{
	return getGenotype(Genotype(c0, c1));
}
