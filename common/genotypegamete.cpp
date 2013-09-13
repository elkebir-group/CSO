/*
 * genotypegamete.cpp
 *
 *  Created on: 10-may-2009
 *      Author: M. El-Kebir
 *
 */

#include "genotypegamete.h"

void GenotypeGamete::computeGametes(int nLoci, const DoubleMatrix& RM, double probLowerBound)
{
	_gametes.clear();

	std::vector<int> homozygousLoci;
	std::vector<int> heterozygousLoci;

	int startVal = chromosomeCompare(nLoci, _c0, _c1, homozygousLoci, heterozygousLoci);

	int heterozygousCount = (int) heterozygousLoci.size();
	int heterozygousDim = 1 << heterozygousCount;

	for (int i = 0; i < heterozygousDim; i++)
	{
		int val = startVal;
		for (int j = 0; j < heterozygousCount; j++)
		{
			// i = number to be mapped
			// j = bit-index
			// heterozygousLoci[j] = mapped bit-index
			val |= ((i >> j) & 1) << heterozygousLoci[j];
		}

		Gamete gamete;
		gamete._c = val;
		gamete._prob = computeProb(nLoci, RM, val, homozygousLoci, heterozygousLoci);

		if (gamete._prob >= probLowerBound)
		{
			_gametes.push_back(gamete);
		}
	}
}

void GenotypeGamete::computeGametesCumulative()
{
	int n = _gametes.size();
	_gametesCumulative = GameteVector(n + 1);
	
	_gametesCumulative[0] = g_InvalidGamete;
	
	double sum = 0;
	int i = 1;
	for (GameteVector::const_iterator it = _gametes.begin(); 
		it != _gametes.end(); it++)
	{
		_gametesCumulative[i++] = *it;
		sum += it->_prob;
	}

	for (int i = 0; i < n; i++)
	{
		_gametesCumulative[i+1]._prob /= sum;
		_gametesCumulative[i+1]._prob += _gametesCumulative[i]._prob;
	}
}

void GenotypeGamete::computeGametes(int nLoci, const DoubleMatrix& RM, const Genotype& ideotype,
	double probLowerBound, bool restrictGametes, int maxCrossOver)
{
	_gametes.clear();

	int minimalSubGroupSize = -1;
	if (restrictGametes)
	{
		minimalSubGroupSize = std::max(largestSubGroupSize(nLoci, _c0, ideotype.getC0()), 
			largestSubGroupSize(nLoci, _c1, ideotype.getC0()));
		
		if (!ideotype.isHomozygous())
		{
			minimalSubGroupSize = std::min(minimalSubGroupSize, 
				std::max(largestSubGroupSize(nLoci, _c0, ideotype.getC1()), 
					largestSubGroupSize(nLoci, _c1, ideotype.getC1())));
		}
	}

	std::vector<int> homozygousLoci;
	std::vector<int> heterozygousLoci;

	std::vector<int> crossOverPoints;
	chromosomeCompare(nLoci, _c0, _c1, homozygousLoci, heterozygousLoci);
	for (int i = 0; i < (int) heterozygousLoci.size(); i++)
	{
		int val = heterozygousLoci[i];
		if (val == 0)
		{
			crossOverPoints.push_back(val);
		}
		else if (val == nLoci - 1)
		{
			if (crossOverPoints.size() == 0 || crossOverPoints.back() != val - 1)
				crossOverPoints.push_back(val - 1);
		}
		else
		{
			if (crossOverPoints.size() == 0 || crossOverPoints.back() != val - 1)
				crossOverPoints.push_back(val - 1);
			crossOverPoints.push_back(val);
		}
	}

	int addedGametesCount = 0;

	/* add the two gametes that can be obtained without any cross-over */
	addedGametesCount = addGamete(nLoci, RM, ideotype, probLowerBound, 
		minimalSubGroupSize, homozygousLoci, heterozygousLoci, _c0);
	addedGametesCount += addGamete(nLoci, RM, ideotype, probLowerBound, 
		minimalSubGroupSize, homozygousLoci, heterozygousLoci, _c1);

	if (maxCrossOver >= 1 && addedGametesCount > 0)
	{
		addedGametesCount = computeGametesOneCrossOver(nLoci, RM, ideotype, 
			probLowerBound, minimalSubGroupSize, 
			homozygousLoci, heterozygousLoci, crossOverPoints);
	}
	if (maxCrossOver >= 2 && addedGametesCount > 0)
	{
		addedGametesCount = computeGametesTwoCrossOvers(nLoci, RM, ideotype, 
			probLowerBound, minimalSubGroupSize,
			homozygousLoci, heterozygousLoci, crossOverPoints);
	}
	else if (maxCrossOver >= 3 && addedGametesCount > 0)
	{
		for (int i = 3; i <= std::min(maxCrossOver, (int) crossOverPoints.size()) && addedGametesCount > 0; i++)
		{
			addedGametesCount = 0;
			computeGametes(nLoci, RM, ideotype, probLowerBound, minimalSubGroupSize, 
				i, homozygousLoci, heterozygousLoci, crossOverPoints, std::vector<int>(), addedGametesCount);
		}
	}

	/* remove duplicates */
	std::sort(_gametes.begin(), _gametes.end(), gamete_lt);
	_gametes.erase(std::unique(_gametes.begin(), _gametes.end(), gamete_eq), _gametes.end());
}

/* This is the generalization of computeGametesOneCrossOver and computeGametesTwoCrossOvers */
void GenotypeGamete::computeGametes(int nLoci, const DoubleMatrix& RM, const Genotype& ideotype,
	double probLowerBound, int minimalSubGroupSize, int crossOverCount, 
	const std::vector<int>& homozygousLoci, const std::vector<int>& heterozygousLoci,
	const std::vector<int>& availableCrossOverPoints, const std::vector<int>& fixedCrossOverPoints,
	int& addedGametesCount)
{
	if (crossOverCount == 0)
	{
		std::vector<int> cpyFixedCrossOverPoints;
		cpyFixedCrossOverPoints.push_back(-1);
		cpyFixedCrossOverPoints.insert(cpyFixedCrossOverPoints.end(), fixedCrossOverPoints.begin(), fixedCrossOverPoints.end());
		cpyFixedCrossOverPoints.push_back(nLoci - 1);

		int mask = 0;
		int n = cpyFixedCrossOverPoints.size() - 1;
		for (int i = 0; i < n; i += 2)
		{
			mask |= GENERATE_AND_MASK(cpyFixedCrossOverPoints[i] + 1, 
				cpyFixedCrossOverPoints[i+1] - cpyFixedCrossOverPoints[i]);
		}

		int complementMask = (~mask) & (GENERATE_AND_MASK(0, nLoci));

		int gamete1 = (_c0 & mask) | (_c1 & complementMask);
		int gamete2 = (_c1 & mask) | (_c0 & complementMask);

		addedGametesCount += addGamete(nLoci, RM, ideotype, probLowerBound, 
			minimalSubGroupSize, homozygousLoci, heterozygousLoci, gamete1) ? 1 : 0;

		addedGametesCount += addGamete(nLoci, RM, ideotype, probLowerBound, 
			minimalSubGroupSize, homozygousLoci, heterozygousLoci, gamete2) ? 1 : 0;
	}
	else
	{
		for (std::vector<int>::const_iterator it = availableCrossOverPoints.begin();
			it != availableCrossOverPoints.end(); it++)
		{
			std::vector<int> newFixedCrossOverPoints = fixedCrossOverPoints;
			newFixedCrossOverPoints.push_back(*it);

			std::vector<int> newAvailableCrossOverPoints;
			newAvailableCrossOverPoints.insert(newAvailableCrossOverPoints.end(), it + 1, availableCrossOverPoints.end());

			computeGametes(nLoci, RM, ideotype,	probLowerBound, minimalSubGroupSize, crossOverCount - 1, 
				homozygousLoci, heterozygousLoci,
				newAvailableCrossOverPoints, newFixedCrossOverPoints, addedGametesCount);
		}
	}
}

int GenotypeGamete::computeGametesOneCrossOver(int nLoci, const DoubleMatrix& RM, const Genotype& ideotype,
	double probLowerBound, int minimalSubGroupSize, const std::vector<int>& homozygousLoci, 
	const std::vector<int>& heterozygousLoci, const std::vector<int>& availableCrossOverPoints)
{
	int nGametesAdded = 0;

	for (int i = 0; i < (int) availableCrossOverPoints.size(); i++)
	{
		int mask = GENERATE_AND_MASK(0, availableCrossOverPoints[i] + 1);
		int compMask = (~mask) & (GENERATE_AND_MASK(0, nLoci));

		int gamete0 = (_c0 & mask) | (_c1 & compMask);
		int gamete1 = (_c1 & mask) | (_c0 & compMask);

		nGametesAdded += addGamete(nLoci, RM, ideotype, probLowerBound, minimalSubGroupSize,
			homozygousLoci, heterozygousLoci, gamete0) ? 1 : 0;
		nGametesAdded += addGamete(nLoci, RM, ideotype, probLowerBound, minimalSubGroupSize, 
			homozygousLoci, heterozygousLoci, gamete1) ? 1 : 0;
	}

	return nGametesAdded;
}

int GenotypeGamete::computeGametesTwoCrossOvers(int nLoci, const DoubleMatrix& RM, const Genotype& ideotype,
	double probLowerBound, int minimalSubGroupSize, const std::vector<int>& homozygousLoci, 
	const std::vector<int>& heterozygousLoci, const std::vector<int>& availableCrossOverPoints)
{
	int nGametesAdded = 0;

	for (int i = 0; i < (int) availableCrossOverPoints.size(); i++)
	{
		for (int j = i + 1; j < (int) availableCrossOverPoints.size(); j++)
		{
			int mask = GENERATE_AND_MASK(0, availableCrossOverPoints[i] + 1) | 
				GENERATE_AND_MASK(availableCrossOverPoints[j] + 1, nLoci - (availableCrossOverPoints[j] + 1));
			int compMask = (~mask) & (GENERATE_AND_MASK(0, nLoci));

			int gamete0 = (_c0 & mask) | (_c1 & compMask);
			int gamete1 = (_c1 & mask) | (_c0 & compMask);

			nGametesAdded += addGamete(nLoci, RM, ideotype, probLowerBound, minimalSubGroupSize,
				homozygousLoci, heterozygousLoci, gamete0) ? 1 : 0;
			nGametesAdded += addGamete(nLoci, RM, ideotype, probLowerBound, minimalSubGroupSize, 
				homozygousLoci, heterozygousLoci, gamete1) ? 1 : 0;
		}
	}

	return nGametesAdded;
}
