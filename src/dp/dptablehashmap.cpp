/*
 * dptablehashmap.h
 *
 *  Created on: 09-mar-2009
 *      Author: s030858
 */

#include "dptablehashmap.h"
#include "data.h"

DpTableHashMap::DpTableHashMap(int nLoci)
	: DpTable(nLoci)
	, _M(1 << nLoci)
{
}

DpTableHashMap::~DpTableHashMap()
{
}

bool DpTableHashMap::isPresent(const Genotype& genotype) const
{
	int c0 = genotype.getC0();
	int c1 = genotype.getC1();

	DpMatrix::const_iterator it = _M.find(c0);
	if (it == _M.end()) return false;
	return it->second.find(c1) != it->second.end();
}

const DpItem& DpTableHashMap::getItem(const Genotype& genotype) const
{
	static const int nLoci = Data::getInstance()->getNumberOfLoci();
	static const DoubleMatrix& RM = Data::getInstance()->getRM();
	static const double probLowerBound = Data::getInstance()->getProbLowerBound();

	int c0 = genotype.getC0();
	int c1 = genotype.getC1();

	DpMatrix::iterator it = _M.find(c0);
	if (it == _M.end())
	{
		DpItem item(genotype, NULL, NULL, INT_MAX, INT_MAX, INT_MAX, INT_MAX, INT_MAX, GenotypeSet(), false);
		item.computeGametes(nLoci, RM, probLowerBound);

		_useCount++;
		
		return _M[c0].insert(std::make_pair(c1, item)).first->second;
	}
	else 
	{
		std::tr1::unordered_map<int, DpItem>::iterator it2 = it->second.find(c1);
		if (it2 == it->second.end())
		{
			DpItem item(genotype, NULL, NULL, INT_MAX, INT_MAX, INT_MAX, INT_MAX, INT_MAX, GenotypeSet(), false);
			item.computeGametes(nLoci, RM, probLowerBound);

			_useCount++;

			return it->second.insert(std::make_pair(item.getC1(), item)).first->second;
		}
		else
		{
			return it2->second;
		}
	}
}

DpItem& DpTableHashMap::getItem(const Genotype& genotype)
{
	static const int nLoci = Data::getInstance()->getNumberOfLoci();
	static const DoubleMatrix& RM = Data::getInstance()->getRM();
	static const double probLowerBound = Data::getInstance()->getProbLowerBound();

	int c0 = genotype.getC0();
	int c1 = genotype.getC1();

	DpMatrix::iterator it = _M.find(c0);
	if (it == _M.end())
	{
		DpItem item(genotype, NULL, NULL, INT_MAX, INT_MAX, INT_MAX, INT_MAX, INT_MAX, GenotypeSet(), false);
		item.computeGametes(nLoci, RM, probLowerBound);

		_useCount++;
		
		return _M[c0].insert(std::make_pair(c1, item)).first->second;
	}
	else 
	{
		std::tr1::unordered_map<int, DpItem>::iterator it2 = it->second.find(c1);
		if (it2 == it->second.end())
		{
			DpItem item(genotype, NULL, NULL, INT_MAX, INT_MAX, INT_MAX, INT_MAX, INT_MAX, GenotypeSet(), false);
			item.computeGametes(nLoci, RM, probLowerBound);

			_useCount++;

			return it->second.insert(std::make_pair(item.getC1(), item)).first->second;
		}
		else
		{
			return it2->second;
		}
	}
}

DpItem& DpTableHashMap::updateItem(const DpItem& item)
{
	DpItem& itemToUpdate = getItem(item);
	itemToUpdate.update(item);
	return itemToUpdate;
}
