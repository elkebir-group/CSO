/*
 * coveranalysis.cpp
 *
 *  Created on: 04-may-2009
 *      Author: M. El-Kebir
 */

#include "coveranalysis.h"

CoverAnalysis::CoverAnalysis(const Data* pData)
{
	const int nLoci = pData->getNumberOfLoci();
	const Genotype& ideotype = pData->getIdeotype();
	const GenotypeSet& parents = pData->getParents();
	const DoubleMatrix& RM = pData->getRM();
	const double gamma = pData->getGamma();

	analyse(nLoci, parents, ideotype.getC0(), RM, gamma);
	
	if (!ideotype.isHomozygous())
		analyse(nLoci, parents, ideotype.getC1(), RM, gamma);
}

CoverAnalysis::CoverAnalysis(const int nLoci, const GenotypeSet& parents, 
	const Genotype& ideotype, const DoubleMatrix& RM, const double gamma)
{
	analyse(nLoci, parents, ideotype.getC0(), RM, gamma);
	
	if (!ideotype.isHomozygous())
		analyse(nLoci, parents, ideotype.getC1(), RM, gamma);
}

void CoverAnalysis::analyse(const int nLoci, const GenotypeSet& parents, const int targetChromosome, 
	const DoubleMatrix& RM, const double gamma)
{
	GenotypePopVector& genotypePopVector = _genotypePopVectorMap[targetChromosome];
	for (int i = 0; i < nLoci; i++)
	{
		GenotypePopEntry entry = { Genotype(0, 0), ULONG_MAX, 0 };
		genotypePopVector.push_back(entry);
	}

	for (int i = 0; i < nLoci; i++)
	{
		for (GenotypeSet::const_iterator it = parents.begin(); it != parents.end(); it++)
		{
            if ((*it)(0, i) == GET_BIT(nLoci, targetChromosome, i))
			{
				genotypePopVector[i]._count++;
				
				unsigned long pop = probToPop(it->computeProb(nLoci, RM, it->getC0()), gamma);
				if (pop <= genotypePopVector[i]._pop)
				{
					genotypePopVector[i]._genotype = *it;
					genotypePopVector[i]._pop = pop;
				}
			}
            else if ((*it)(1, i) == GET_BIT(nLoci, targetChromosome, i))
			{
				genotypePopVector[i]._count++;
				
				unsigned long pop = probToPop(it->computeProb(nLoci, RM, it->getC1()), gamma);
				if (pop <= genotypePopVector[i]._pop)
				{
					genotypePopVector[i]._genotype = *it;
					genotypePopVector[i]._pop = pop;
				}
			}
		}
	}
}

unsigned long CoverAnalysis::getPopMax() const
{
	unsigned long max = 0;
	for (GenotypePopVectorMap::const_iterator it = _genotypePopVectorMap.begin(); 
		it != _genotypePopVectorMap.end(); it++)
	{
		const GenotypePopVector& popVector = it->second;
		for (GenotypePopVector::const_iterator it2 = popVector.begin(); it2 != popVector.end(); it2++)
		{
			if (max < it2->_pop)
			{
				max = it2->_pop;
			}
		}
	}

	return max;
}

const GenotypePopVector& CoverAnalysis::getGenotypePopVector(const int targetChromosome) const
{
	assert(_genotypePopVectorMap.find(targetChromosome) != _genotypePopVectorMap.end());

	return _genotypePopVectorMap.find(targetChromosome)->second;
}


GenotypeSet CoverAnalysis::getMinimalCover() const
{
	GenotypeSet cover;
	for (GenotypePopVectorMap::const_iterator it = _genotypePopVectorMap.begin(); 
		it != _genotypePopVectorMap.end(); it++)
	{
		const GenotypePopVector& popVector = it->second;
		for (GenotypePopVector::const_iterator it2 = popVector.begin(); it2 != popVector.end(); it2++)
		{
			if (it2->_count == 1)
			{
				cover.insert(it2->_genotype);
			}
		}
	}

	return cover;
}
