/*
 * cso.cpp
 *
 *  Created on: 18-feb-2009
 *      Author: M. El-Kebir
 */

/* Testing
 
#include "cso.h"

const Gamete g_InvalidGamete = {INVALID_GAMETE, 0};

void toBitstring(int val, int n, char* buf)
{
	/*
	 * Example:
	 * '1000' = 7
	 * '0001' = 1
	 */
	for (int i = 0; i < n; i++)
	{
		buf[n - i - 1] = (val & 1) ? '1' : '0';
		val >>= 1;
	}

	buf[n] = '\0';
}

int fromBitstring(int n, const char* buf)
{
	/*
	 * Example:
	 * '1000' = 8
	 * '0001' = 1
	 */
	int val = 0;
	for (int i = 0; i < n; i++)
	{
		if (buf[n - i - 1] == '1')
		{
			val |= 1 << i;
		}
	}
	return val;
}

int numberOfDifferences(int n, int val1, int val2)
{
	int res = 0;

	for (int i = 0; i < n; i++)
	{
		if ((val1 & 1) != (val2 & 1)) res++;

		val1 >>= 1;
		val2 >>= 1;
	}

	return res;
}

/*
 * Index 0 = LSB
 */
int chromosomeCompare(int n, int val1, int val2,
		std::vector<int>& homyzygousLoci, std::vector<int>& heterozygousLoci)
{
	homyzygousLoci.clear();
	heterozygousLoci.clear();

	int base = 0;
	for (int i = 0; i < n; i++)
	{
		if ((val1 & 1) == (val2 & 1))	// homozygous
		{
			base |= (val1 & 1) << i;
			homyzygousLoci.push_back(i);
		}
		else				// heterozygous
			heterozygousLoci.push_back(i);

		val1 >>= 1;
		val2 >>= 1;
	}

	return base;
}

int randInt(int n)
{
	return (int) (n * (rand() / ((double)RAND_MAX + 1)));
}

int randChromosome(int nLoci)
{
	int res = 0;
	for (int i = 0; i < nLoci; i++)
	{
		if (randDouble() < 0.5)
			res |= 1 << i;
	}

	return res;
}

double randDouble()
{
	return rand() / ((double)RAND_MAX + 1);
}

unsigned long probToPop(double p, double gamma)
{
	if ((1 - p) <= (2 * DBL_EPSILON)) return 1;
        return (unsigned long) ceil(log(1 - gamma) / log(1 - p));
}

int largestSubGroupSize(int nLoci, int val, int target)
{
	int largestGroupSize = 0;
	int currentGroupSize = 0;

	for (int i = 0; i < nLoci; i++)
	{
		if (GET_BIT(val, i) == GET_BIT(target, i))
		{
			currentGroupSize++;
			if (currentGroupSize > largestGroupSize)
			{
				largestGroupSize = currentGroupSize;
			}
		}
		else
		{
			currentGroupSize = 0;
		}
	}

	return largestGroupSize;
}

bool gamete_eq(const Gamete& gamete1, const Gamete& gamete2)
{
	return gamete1._c == gamete2._c; 
}

bool gamete_lt(const Gamete& gamete1, const Gamete& gamete2)
{
	return gamete1._c < gamete2._c; 
}
