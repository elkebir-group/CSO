/*
 * genotypegamete.h
 *
 *  Created on: 13-apr-2009
 *      Author: s030858
 *
 * Remark: set GameteVector to GameteList for improved perfomance...
 */

#ifndef GENOTYPEGAMETE_H_
#define GENOTYPEGAMETE_H_

#include "cso.h"
#include "genotype.h"
#include <iostream>
#include <math.h>

class GenotypeGamete : public Genotype
{
private:
	GameteVector _gametes;
	GameteVector _gametesCumulative;

protected:
	void getGameteList(int nLoci, const DoubleMatrix& RM, double probLowerBound, GameteList& gameteList) const;

public:
	GenotypeGamete(const Genotype& genotype);
	GenotypeGamete(int c0, int c1);
	virtual ~GenotypeGamete();
	void computeGametes(int nLoci, const DoubleMatrix& RM, double probLowerBound);
	void computeGametesCumulative();
	const GameteVector& getGametes() const;
	Gamete getGameteMostAlike(int nLoci, int target) const;
	Gamete getGameteMostAlikeProportionate(int nLoci, int target) const;
	Gamete getGameteUniformlyAtRandom() const;
	Gamete getGameteProportionate() const;
};

inline GenotypeGamete::GenotypeGamete(const Genotype& genotype)
	: Genotype(genotype)
	, _gametes()
	, _gametesCumulative()
{
}

inline GenotypeGamete::GenotypeGamete(int c0, int c1)
	: Genotype(c0, c1)
	, _gametes()
	, _gametesCumulative()
{
}

inline GenotypeGamete::~GenotypeGamete()
{
}

inline void GenotypeGamete::computeGametes(int nLoci, const DoubleMatrix& RM, double probLowerBound)
{
	_gametes.clear();

	std::vector<int> homozygousLoci;
	std::vector<int> heterozygousLoci;

	int startVal = chromosomeCompare(nLoci, _c0, _c1, homozygousLoci, heterozygousLoci);

	int heterozygousCount = (int) heterozygousLoci.size();
	int heterozygousDim = 1 << heterozygousCount;

	//double max = DBL_MAX;

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

			/*if (max > gamete._prob)
			{
				max = gamete._prob;

				res.clear();
				res.push_back(gamete);
			}*/
		}
	}
}

inline void GenotypeGamete::computeGametesCumulative()
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

inline const GameteVector& GenotypeGamete::getGametes() const
{
	return _gametes;
}

inline Gamete GenotypeGamete::getGameteMostAlike(int nLoci, int target) const
{
	int min = INT_MAX;
	Gamete res = g_InvalidGamete;
	
	for (GameteVector::const_iterator it = _gametes.begin(); 
		it != _gametes.end(); it++)
	{
		int diff = numberOfDifferences(nLoci, target, it->_c);
		if (diff < min)
			res = *it;
	}

	return res;
}

inline Gamete GenotypeGamete::getGameteMostAlikeProportionate(int nLoci, int target) const
{
	int n = _gametes.size();
	GameteVector gametesCumulative(n + 1);
	
	gametesCumulative[0] = g_InvalidGamete;
	
	double sum = 0;
	int i = 1;
	for (GameteVector::const_iterator it = _gametes.begin(); 
		it != _gametes.end(); it++)
	{
		gametesCumulative[i] = *it;
		//gametesCumulative[i]._prob /= (1 + numberOfDifferences(nLoci, target, it->_c));
		gametesCumulative[i]._prob = 1.0 / (numberOfDifferences(nLoci, target, it->_c) + 1);
		sum += gametesCumulative[i++]._prob;
	}

	for (int i = 0; i < n; i++)
	{
		gametesCumulative[i+1]._prob /= sum;
		gametesCumulative[i+1]._prob += gametesCumulative[i]._prob;
	}

	double p = randDouble();

	/* binary search, invariants: 
	 * - f[x] \leq A < f[y] 
	 * - 0 \leq x < y \leq N
	 */
	int x = 0, y = (int) gametesCumulative.size();

	if (y == 1)
		return g_InvalidGamete;

	while (x + 1 != y)
	{
		// can overflow: x + ((y - x) / 2) is safer
		int h = (x + y) / 2;
		
		if (gametesCumulative[h]._prob < p)
			x = h;
		else
			y = h;
	}

	return _gametes[x];
}

inline Gamete GenotypeGamete::getGameteUniformlyAtRandom() const
{
	Gamete res = g_InvalidGamete;
	
	int nGametes = (int) _gametes.size();
	if (nGametes)
	{
		res = _gametes[randInt((int) _gametes.size())];
	}

	return res;
}

inline Gamete GenotypeGamete::getGameteProportionate() const
{
	double p = randDouble();

	/* binary search, invariants: 
	 * - f[x] \leq A < f[y] 
	 * - 0 \leq x < y \leq N
	 */
	int x = 0, y = (int) _gametesCumulative.size();

	if (y == 1)
		return g_InvalidGamete;

	while (x + 1 != y)
	{
		// can overflow: x + ((y - x) / 2) is safer
		int h = (x + y) / 2;
		
		if (_gametesCumulative[h]._prob < p)
			x = h;
		else
			y = h;
	}

	//Gamete res = g_InvalidGamete;
	//if (x >= 0)
	return _gametes[x];
}

#endif