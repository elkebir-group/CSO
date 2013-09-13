/*
 * genotypetable.h
 *
 *  Created on: 20-mar-2009
 *      Author: M. El-Kebir
 */

#ifndef GENOTYPETABLE_H_
#define GENOTYPETABLE_H_

#include "csoga.h"
#include "genotypegamete.h"
#include "genotype.h"

class GenotypeTable
{
protected:
	const int _nLoci;
	const int _nGenotypes;
	const DoubleMatrix _RM;
	const double _probLowerBound;
	const Genotype _ideotype;
	const unsigned int _limit;
	const int _crossOverCount;
	const bool _restrictGametes;
	mutable GenotypeList _genotypeList;
	static GenotypeTable* _pGenotypeTable;

protected:
	GenotypeTable(unsigned int limit, int crossOverCount, bool restrictGametes);

public:
	virtual ~GenotypeTable();
	virtual const GenotypeGamete& getGenotype(const Genotype& genotype) const = 0;
	virtual const GenotypeGamete& getGenotype(int c0, int c1) const = 0;
	unsigned int getUseCount() const;
	int getCrossOverCount() const;
	bool getRestrictGametes() const;
	static const GenotypeTable* getInstance();
};

inline const GenotypeTable* GenotypeTable::getInstance()
{
	return _pGenotypeTable;
}

inline unsigned int GenotypeTable::getUseCount() const
{
	return _genotypeList.size();
}

inline int GenotypeTable::getCrossOverCount() const
{
	return _crossOverCount;
}

inline bool GenotypeTable::getRestrictGametes() const
{
	return _restrictGametes;
}

#endif
