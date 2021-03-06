/*
 * cso.cpp
 *
 *  Created on: 18-feb-2009
 *      Author: s030858
 */

#include "cso.h"

const Gamete g_InvalidGamete = {INVALID_GAMETE, 0};

void toBitstring(int val, int n, char* buf)
{
	/*
	 * Example:
	 * '1000' = 8
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
      val |= (1 << i);
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
    if (GET_BIT(n, val1, i) == GET_BIT(n, val2, i))	// homozygous
		{
      base |= GET_BIT(n, val1, i) << (n - i - 1);
      homyzygousLoci.push_back(i);
    }
		else				// heterozygous
      heterozygousLoci.push_back(i);
	}

  // base is the case where all bits at heterozygous loci are set to 0
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
