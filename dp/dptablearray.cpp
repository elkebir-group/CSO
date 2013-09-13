/*
 * dptablearray.h
 *
 *  Created on: 09-mar-2009
 *      Author: M. El-Kebir
 */

#include "dptablearray.h"
#include "data.h"
#include <string.h>

DpTableArray::DpTableArray(int crossOverCount, bool restrictGametes)
	: DpTable(crossOverCount, restrictGametes)
	, _M(new DpItem*[_nGenotypes])
{
	memset(_M, NULL, sizeof(DpItem*) * _nGenotypes);
}

DpTableArray::~DpTableArray()
{
	for (int i = 0; i < _nGenotypes; i++)
	{
		delete _M[i];
	}

	delete[] _M;
}

bool DpTableArray::isPresent(const Genotype& genotype) const
{
	assert(getIndex(genotype) < _nGenotypes);
	return _M[getIndex(genotype)] != NULL;
}

const DpItem& DpTableArray::getItem(const Genotype& genotype) const
{
	int index = getIndex(genotype);
	
	if (_M[index])
	{
		return *_M[getIndex(genotype)];
	}
	else
	{
		_useCount++;

		_M[index] = new DpItem(genotype, NULL, NULL, INT_MAX, INT_MAX, INT_MAX, INT_MAX, GenotypeSet(), false);

		if (_crossOverCount > 0) 
			_M[index]->computeGametes(_nLoci, _RM, _ideotype, 
				_probLowerBound, _restrictGametes, _crossOverCount);
		else
			_M[index]->computeGametes(_nLoci, _RM, _probLowerBound);

		return *_M[index];
	}
}

DpItem& DpTableArray::getItem(const Genotype& genotype)
{
	int index = getIndex(genotype);
	
	if (_M[index])
	{
		return *_M[getIndex(genotype)];
	}
	else
	{
		_useCount++;

		_M[index] = new DpItem(genotype, NULL, NULL, INT_MAX, INT_MAX, INT_MAX, INT_MAX, GenotypeSet(), false);

		if (_crossOverCount > 0) 
			_M[index]->computeGametes(_nLoci, _RM, _ideotype, 
				_probLowerBound, _restrictGametes, _crossOverCount);
		else
			_M[index]->computeGametes(_nLoci, _RM, _probLowerBound);

		return *_M[index];
	}
}

DpItem& DpTableArray::updateItem(const DpItem& item)
{
	DpItem& itemToUpdate = getItem(item);
	itemToUpdate.update(item);
	return itemToUpdate;
}
