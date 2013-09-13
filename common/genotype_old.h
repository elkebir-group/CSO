/*
 * genotype.h
 *
 *  Created on: 19-feb-2009
 *      Author: M. El-Kebir
 *
 * Remark: set GameteVector to GameteList for improved perfomance...
 */

#ifndef GENOTYPE_H_
#define GENOTYPE_H_

#include "cso.h"
#include <iostream>
#include <math.h>

class Genotype
{
protected:
	int _c0;
	int _c1;
	double computeProb(int nLoci, const DoubleMatrix& RM, int gamete,
			const std::vector<int>& homozygousLoci, const std::vector<int>& heterozygousLoci) const;
	static int getAllele(int locus, int chromosome);

public:
	Genotype(int c0, int c1);
	virtual ~Genotype();
	int getC0() const;
	int getC1() const;
	int compare(const Genotype& genotype) const;
	bool operator <(const Genotype& genotype) const;
	bool operator <=(const Genotype& genotype) const;
	bool operator >(const Genotype& genotype) const;
	bool operator >=(const Genotype& genotype) const;
	bool operator ==(const Genotype& genotype) const;
	bool operator !=(const Genotype& genotype) const;
	int operator ()(const int i, const int j) const;
	double computeProb(int nLoci, const DoubleMatrix& RM, int c) const;
	void printGenotype(int nLoci, bool newline = true, std::ostream& out = std::cout, const char* separator = "/") const;
	unsigned long computePop(int nLoci, const DoubleMatrix& RM, double gamma, const Genotype& D, const Genotype& E) const;
	bool isHomozygous() const;
	int getNumberOfHomozygousLoci(int nLoci) const;
	LinkageType getLinkage(int locus1, int locus2, int targetChromosome) const;
	int getMaskAtMostOneCrossOver(int nLoci, int target) const;
	int getMaskAtMostOneCrossOverGroup(int nLoci, int target) const;
	std::vector<int> getMaskVectorAtMostOneCrossOverGroup(int nLoci, int target) const;
	int getMask(int nLoci, int target) const;
};

inline Genotype::Genotype(int c0, int c1)
	: _c0(c0 < c1 ? c0 : c1)
	, _c1(c0 > c1 ? c0 : c1)
{
}

inline Genotype::~Genotype()
{
}

inline int Genotype::getC0() const
{
	return _c0;
}

inline int Genotype::getC1() const
{
	return _c1;
}

inline int Genotype::compare(const Genotype& genotype) const
{
	if (_c0 < genotype._c0)
		return -1;
	else if (_c0 > genotype._c0)
		return 1;
	else if (_c1 < genotype._c1)	// _c0 == genotype._c0
		return -1;
	else if (_c1 > genotype._c1)
		return 1;
	else
		return 0;					// _c1 == genotype._c1
}

inline bool Genotype::operator <(const Genotype& genotype) const
{
	return compare(genotype) < 0;
}

inline bool Genotype::operator <=(const Genotype& genotype) const
{
	return compare(genotype) <= 0;
}

inline bool Genotype::operator >(const Genotype& genotype) const
{
	return compare(genotype) > 0;
}

inline bool Genotype::operator >=(const Genotype& genotype) const
{
	return compare(genotype) >= 0;
}

inline bool Genotype::operator ==(const Genotype& genotype) const
{
	return compare(genotype) == 0;
}

inline bool Genotype::operator !=(const Genotype& genotype) const
{
	return compare(genotype) != 0;
}

inline double Genotype::computeProb(int nLoci, const DoubleMatrix& RM, int gamete,
		const std::vector<int>& homozygousLoci, const std::vector<int>& heterozygousLoci) const
{
	for (std::vector<int>::const_iterator it = homozygousLoci.begin();
		it != homozygousLoci.end(); it++)
	{
		assert((*this)(0, *it) == (*this)(1, *it));

		if (((gamete >> *it) & 1) != ((_c0 >> *it) & 1)) return 0;
	}

	int heterozygousLociCount = (int) heterozygousLoci.size();
	if (heterozygousLociCount == 0)
		return 1;

	double res = 0.5;
	for (int i = 0; i < heterozygousLociCount - 1; i++)
	{
		// val_chromosome_locus
		int val_0_0 = (*this)(0, heterozygousLoci[i]);
		int val_0_1 = (*this)(0, heterozygousLoci[i+1]);
#ifndef NDEBUG		
		int val_1_0 = (*this)(1, heterozygousLoci[i]);
#endif
		int val_1_1 = (*this)(1, heterozygousLoci[i+1]);

		int gamete_val_0 = (gamete >> heterozygousLoci[i]) & 1;
		int gamete_val_1 = (gamete >> heterozygousLoci[i+1]) & 1;

		if (gamete_val_0 == val_0_0)
		{
			if (gamete_val_1 == val_0_1)
			{
				res *= 1 - RM[heterozygousLoci[i]][heterozygousLoci[i + 1]];
			}
			else
			{
				assert(gamete_val_1 == val_1_1);
				res *= RM[heterozygousLoci[i]][heterozygousLoci[i + 1]];
			}
		}
		else
		{
			assert(gamete_val_0 == val_1_0);

			if (gamete_val_1 == val_1_1)
			{
				res *= 1 - RM[heterozygousLoci[i]][heterozygousLoci[i + 1]];
			}
			else
			{
				assert(gamete_val_1 == val_0_1);
				res *= RM[heterozygousLoci[i]][heterozygousLoci[i + 1]];
			}
		}
	}

	return res;
}

inline int Genotype::getAllele(int locus, int chromosome)
{
	return (chromosome >> locus) & 1;
}

inline int Genotype::operator ()(const int i, const int j) const
{
	assert(i == 0 || i == 1);
	if (i == 0)
	{
		return ((_c0 >> j) & 1);
	}
	else
	{
		return ((_c1 >> j) & 1);
	}
}

inline bool Genotype::isHomozygous() const
{
	return _c0 == _c1;
}

inline int Genotype::getNumberOfHomozygousLoci(int nLoci) const
{
	int res = 0;
	for (int i = 0; i < nLoci; i++)
	{
		if ((*this)(0, i) == (*this)(1, i))
			res++;
	}

	return res;
}

inline LinkageType Genotype::getLinkage(int locus1, int locus2, int targetChromosome) const
{
	bool c0_locus1 = (getAllele(locus1, _c0) == getAllele(locus1, targetChromosome));
	bool c0_locus2 = (getAllele(locus2, _c0) == getAllele(locus2, targetChromosome));
	bool c1_locus1 = (getAllele(locus1, _c1) == getAllele(locus1, targetChromosome));
	bool c1_locus2 = (getAllele(locus2, _c1) == getAllele(locus2, targetChromosome));

	if ((c0_locus1 && c0_locus2) || (c1_locus1 && c1_locus2))
	{
		return Linked;
	}
	else if ((c0_locus1 && c1_locus2) || (c1_locus1 && c0_locus2))
	{
		return WeaklyLinked;
	}
	else
	{
		return Unlinked;
	}
}

#endif /* GENOTYPE_H_ */
