/*
 * coveranalysis.h
 *
 *  Created on: 04-may-2009
 *      Author: M. El-Kebir
 */

#ifndef COVERANALYSIS_H_
#define COVERANALYSIS_H_

#include "cso.h"
#include "data.h"
#include <vector>
#include <map>
#include "genotype.h"

typedef struct
{
  Genotype _genotype;
  unsigned long _pop;
  int _count;
} GenotypePopEntry;

typedef std::vector<GenotypePopEntry> GenotypePopVector;
typedef std::map<int, GenotypePopVector> GenotypePopVectorMap;

class CoverAnalysis
{
private:
  GenotypePopVectorMap _genotypePopVectorMap;
  void analyse(const int nLoci, const GenotypeSet& parents, const int targetChromosome,
    const DoubleMatrix& RM, const double gamma);

public:
  CoverAnalysis(const Data* pData);
  CoverAnalysis(const int nLoci, const GenotypeSet& parents,
    const Genotype& ideotype, const DoubleMatrix& RM, const double gamma);
  unsigned long getPopMax() const;
  const GenotypePopVector& getGenotypePopVector(const int targetChromosome) const;
  GenotypeSet getMinimalCover() const;
};

#endif
