/*
 * genotypetablehashmap.h
 *
 *  Created on: 20-mar-2009
 *      Author: M. El-Kebir
 */

#ifndef GENOTYPETABLEHASHMAP_H_
#define GENOTYPETABLEHASHMAP_H_

#include "csoga.h"
#include "genotypetable.h"

class GenotypeTableHashMap : public GenotypeTable
{
private:
	mutable GenotypeGameteMatrix _A;

public:
	GenotypeTableHashMap(unsigned int limit, int crossOverCount = -1, bool restrictGametes = false);
	~GenotypeTableHashMap();
	const GenotypeGamete& getGenotype(const Genotype& genotype) const;
	const GenotypeGamete& getGenotype(int c0, int c1) const;
};

#endif
