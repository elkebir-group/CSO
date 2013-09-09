/*
 * genotypetablearray.h
 *
 *  Created on: 20-mar-2009
 *      Author: M. El-Kebir
 */

#ifndef GENOTYPETABLEARRAY_H_
#define GENOTYPETABLEARRAY_H_

#include "csoga.h"
#include "genotypetable.h"

class GenotypeTableArray : public GenotypeTable
{
private:
	mutable GenotypeGamete** _A;
	int getIndex(const Genotype& genotype) const;
	int getIndex(int c0, int c1) const;

public:
	GenotypeTableArray(unsigned int limit, int crossOverCount = -1, bool restrictGametes = false);
	~GenotypeTableArray();
	const GenotypeGamete& getGenotype(const Genotype& genotype) const;
	const GenotypeGamete& getGenotype(int c0, int c1) const;
};

inline int GenotypeTableArray::getIndex(const Genotype& genotype) const
{
	int i = genotype.getC0();
	int j = genotype.getC1();

	return (1 << _nLoci) * i - (i + 1) * i / 2 + j;
}

inline int GenotypeTableArray::getIndex(int c0, int c1) const
{
	return (1 << _nLoci) * c0 - (c0 + 1) * c0 / 2 + c1;
}

#endif
