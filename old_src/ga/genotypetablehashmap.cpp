/*
 * genotypetablehashmap.cpp
 *
 *  Created on: 20-mar-2009
 *      Author: M. El-Kebir
 */

#include "genotypetablehashmap.h"
#include "cso.h"

GenotypeTableHashMap::GenotypeTableHashMap(unsigned int limit, int crossOverCount, bool restrictGametes)
	: GenotypeTable(limit, crossOverCount, restrictGametes)
	, _A(1 << _nLoci)
{
}

GenotypeTableHashMap::~GenotypeTableHashMap()
{
}

const GenotypeGamete& GenotypeTableHashMap::getGenotype(const Genotype& genotype) const
{
	int c0 = genotype.getC0();
	int c1 = genotype.getC1();

	return getGenotype(c0, c1);
}

const GenotypeGamete& GenotypeTableHashMap::getGenotype(int c0, int c1) const
{
	GenotypeGameteMatrix::iterator it = _A.find(c0 < c1 ? c0 : c1);
	std::tr1::unordered_map<int, GenotypeGamete>::iterator it2;
	if (it == _A.end())
	{
		Genotype genotype(c0, c1);
		_genotypeList.push_back(genotype);

		GenotypeGamete* pResult = &_A[c0 < c1 ? c0 : c1].insert(std::make_pair(c0 > c1 ? c0 : c1, genotype)).first->second;

		if (_crossOverCount > 0)
			pResult->computeGametes(_nLoci, _RM, _ideotype, 
				_probLowerBound, _restrictGametes, _crossOverCount);
		else
			pResult->computeGametes(_nLoci, _RM, _probLowerBound);
		
		pResult->computeGametesCumulative();

		if (_genotypeList.size() == _limit + 1)
		{
			Genotype oldestGenotype = _genotypeList.front();
			_genotypeList.pop_front();
			_A.find(oldestGenotype.getC0())->second.erase(oldestGenotype.getC1());
		}

		assert(_genotypeList.size() <= _limit);

		return *pResult;
	}
	else if ((it2 = it->second.find(c0 > c1 ? c0 : c1)) == it->second.end())
	{
		Genotype genotype(c0, c1);
		_genotypeList.push_back(genotype);

		GenotypeGamete* pResult = &it->second.insert(std::make_pair(c0 > c1 ? c0 : c1, genotype)).first->second;

		if (_crossOverCount > 0)
			pResult->computeGametes(_nLoci, _RM, _ideotype, 
				_probLowerBound, _restrictGametes, _crossOverCount);
		else
			pResult->computeGametes(_nLoci, _RM, _probLowerBound);
		
		pResult->computeGametesCumulative();
		
		if (_genotypeList.size() == _limit + 1)
		{
			Genotype oldestGenotype = _genotypeList.front();
			_genotypeList.pop_front();
			_A.find(oldestGenotype.getC0())->second.erase(oldestGenotype.getC1());
		}

		assert(_genotypeList.size() <= _limit);

		return *pResult;
	}
	else
	{
		return it2->second;
	}
}
