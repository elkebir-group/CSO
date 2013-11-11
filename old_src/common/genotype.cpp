/*
 * genotype.cpp
 *
 *  Created on: 28-apr-2011
 *      Author: M. El-Kebir
 *
 */

#include "genotype.h"

std::istream& operator >>(std::istream& is, Genotype& genotype)
{
  return is >> genotype._c0 >> genotype._c1;
}

std::ostream& operator <<(std::ostream& os, const Genotype& genotype)
{
  return os << genotype._c0 << " " << genotype._c1;
}

int Genotype::getMask(int nLoci, int target) const
{
	for (int i = 0; i < nLoci; i++)
	{
        int c0_bit = GET_BIT(nLoci, _c0, i);
        int c1_bit = GET_BIT(nLoci, _c1, i);
        int target_bit = GET_BIT(nLoci, target, i);

		if (c0_bit != target_bit && c1_bit != target_bit)
		{
			// flip target_bit
			if (target_bit == 1)
                target &= ~(1 << (nLoci - i - 1));
			else
                target |= (1 << (nLoci - i - 1));
		}
	}

	return target;
}

int Genotype::getMaskAtMostOneCrossOverGroup(int nLoci, int target) const
{
    std::vector<int> largestGroupSize0_0(nLoci); // _c0: left to right
    std::vector<int> largestGroupSize1_0(nLoci); // _c1: right to left
    std::vector<int> largestGroupSize0_1(nLoci); // _c0: right to left
    std::vector<int> largestGroupSize1_1(nLoci); // _c1: left to right

	/* First determine largest group size per chromosome */
    largestGroupSize0_0[0] = GET_BIT(nLoci, _c0, 0) == GET_BIT(nLoci, target, 0) ? 1 : 0;
    largestGroupSize1_1[0] = GET_BIT(nLoci, _c1, 0) == GET_BIT(nLoci, target, 0) ? 1 : 0;

	int largestGroup0 = largestGroupSize0_0[0];
	int largestGroup1 = largestGroupSize1_1[0];

	for (int i = 1; i < nLoci; i++)
	{
        if (GET_BIT(nLoci, _c0, i) == GET_BIT(nLoci, target, i))
		{
			largestGroupSize0_0[i] += largestGroupSize0_0[i - 1] + 1;

			if (largestGroupSize0_0[i] > largestGroup0)
			{
				largestGroup0 = largestGroupSize0_0[i];
			}
		}
		else
		{
			largestGroupSize0_0[i] = 0;
		}

        if (GET_BIT(nLoci, _c1, i) == GET_BIT(nLoci, target, i))
		{
			largestGroupSize1_1[i] += largestGroupSize1_1[i - 1] + 1;

			if (largestGroupSize1_1[i] > largestGroup1)
			{
				largestGroup1 = largestGroupSize1_1[i];
			}
		}
		else
		{
			largestGroupSize1_1[i] = 0;
		}
	}

    largestGroupSize0_1[nLoci - 1] = GET_BIT(nLoci, _c0, nLoci - 1) == GET_BIT(nLoci, target, nLoci - 1) ? 1 : 0;
    largestGroupSize1_0[nLoci - 1] = GET_BIT(nLoci, _c1, nLoci - 1) == GET_BIT(nLoci, target, nLoci - 1) ? 1 : 0;
	for (int i = nLoci - 2; i >= 0; i--)
	{
        if (GET_BIT(nLoci, _c1, i) == GET_BIT(nLoci, target, i))
		{
			largestGroupSize1_0[i] += largestGroupSize1_0[i + 1] + 1;
		}
		else
		{
			largestGroupSize1_0[i] = 0;
		}

        if (GET_BIT(nLoci, _c0, i) == GET_BIT(nLoci, target, i))
		{
			largestGroupSize0_1[i] += largestGroupSize0_1[i + 1] + 1;
		}
		else
		{
			largestGroupSize0_1[i] = 0;
		}
	}

	/* Determine best crossover point */

	int bestCrossoverPoint = -1;
	int largestGroup = std::max(largestGroup0, largestGroup1);
	bool rightToLeftC0 = true;

	for (int i = 0; i < nLoci - 1; i++)
	{
		if ((largestGroupSize0_0[i] + largestGroupSize1_0[i + 1]) > largestGroup)
		{
			largestGroup = largestGroupSize0_0[i] + largestGroupSize1_0[i + 1];
			bestCrossoverPoint = i;
			rightToLeftC0 = true;
		}
		
		if ((largestGroupSize1_1[i] + largestGroupSize0_1[i + 1]) > largestGroup)
		{
			largestGroup = largestGroupSize1_1[i] + largestGroupSize0_1[i + 1];
			bestCrossoverPoint = i;
			rightToLeftC0 = false;
		}
	}

	/* Generate the mask */

	if (bestCrossoverPoint == -1 && largestGroup0 > largestGroup1)
	{
		return _c0;
	}
	else if (bestCrossoverPoint == -1) // largestGroup0 <= largestGroup1
	{
		return _c1;
	}
	else
	{
		int val1 = 0, val2 = 0;
		for (int i = 0; i < nLoci; i++)
		{
			if (i <= bestCrossoverPoint)
			{
				val1 |= (_c0 & (1 << i));
				val2 |= (_c1 & (1 << i));
			}
			else
			{
				val1 |= (_c1 & (1 << i));
				val2 |= (_c0 & (1 << i));
			}
		}

		return rightToLeftC0 ? val1 : val2;
	}
}

int Genotype::getMaskAtMostOneCrossOver(int nLoci, int target) const
{
	/* first determine the best crossover point; there are nLoci - 1 crossover points */
	std::vector<int> correctAlleleCount0(nLoci - 1);
	std::vector<int> correctAlleleCount1(nLoci - 1);
	
	correctAlleleCount0[0] = (target & 1) == (_c0 & 1) ? 1 : 0;
	for (int i = 1; i < nLoci - 1; i++)
	{
		correctAlleleCount0[i] += correctAlleleCount0[i - 1] + (((target >> i) & 1) == ((_c0 >> i) & 1) ? 1 : 0);
	}

	correctAlleleCount1[nLoci - 2] += ((target >> (nLoci - 1)) & 1) == ((_c1 >> (nLoci - 1)) & 1) ? 1 : 0;
	for (int i = nLoci - 3; i >= 0; i--)
	{
		correctAlleleCount1[i] += correctAlleleCount1[i + 1] + (((target >> (i + 1)) & 1) == ((_c1 >> (i + 1)) & 1) ? 1 : 0);
	}

	std::vector<int> correctAlleleCount(nLoci - 1);
	for (int i = 0; i < nLoci - 1; i++)
	{
		correctAlleleCount[i] = correctAlleleCount0[i] + correctAlleleCount1[i];
	}

	int bestCrossOverPoint = -1;
	int bestCorrectAlleleCount = 0;
	for (int i = 0; i < nLoci - 1; i++)
	{
		if (correctAlleleCount[i] > bestCorrectAlleleCount)
		{
			bestCrossOverPoint = i;
			bestCorrectAlleleCount = correctAlleleCount[i];
		}
	}

	/* now generate the mask */
	if (bestCorrectAlleleCount <= (nLoci - numberOfDifferences(nLoci, _c0, target)))
	{
		return _c0;
	}
	else if (bestCorrectAlleleCount <= (nLoci - numberOfDifferences(nLoci, _c1, target)))
	{
		return _c1;
	}
	else
	{
		int val1 = 0;
		int val2 = 0;
		for (int i = 0; i < nLoci; i++)
		{
			if (i <= bestCrossOverPoint)
			{
				val1 |= (_c0 & (1 << i));
				val2 |= (_c1 & (1 << i));
			}
			else
			{
				val1 |= (_c1 & (1 << i));
				val2 |= (_c0 & (1 << i));
			}
		}

		if (numberOfDifferences(nLoci, val1, target) > numberOfDifferences(nLoci, val2, target))
		{
			return val2;
		}
		else
		{
			return val1;
		}
	}
}
std::vector<int> Genotype::getMaskVectorAtMostOneCrossOverGroup(int nLoci, int target) const
{
	std::vector<int> largestGroupSize0_0(nLoci); // _c0: right to left
	std::vector<int> largestGroupSize1_0(nLoci); // _c1: left to right
	std::vector<int> largestGroupSize0_1(nLoci); // _c0: left to right
	std::vector<int> largestGroupSize1_1(nLoci); // _c1: right to left

	/* First determine largest group size per chromosome */
    largestGroupSize0_0[0] = GET_BIT(nLoci, _c0, 0) == GET_BIT(nLoci, target, 0) ? 1 : 0;
    largestGroupSize1_1[0] = GET_BIT(nLoci, _c1, 0) == GET_BIT(nLoci, target, 0) ? 1 : 0;

	int largestGroup0 = largestGroupSize0_0[0];
	int largestGroup1 = largestGroupSize1_1[0];

	for (int i = 1; i < nLoci; i++)
	{
        if (GET_BIT(nLoci,_c0, i) == GET_BIT(nLoci, target, i))
		{
			largestGroupSize0_0[i] += largestGroupSize0_0[i - 1] + 1;

			if (largestGroupSize0_0[i] > largestGroup0)
			{
				largestGroup0 = largestGroupSize0_0[i];
			}
		}
		else
		{
			largestGroupSize0_0[i] = 0;
		}

        if (GET_BIT(nLoci, _c1, i) == GET_BIT(nLoci, target, i))
		{
			largestGroupSize1_1[i] += largestGroupSize1_1[i - 1] + 1;

			if (largestGroupSize1_1[i] > largestGroup1)
			{
				largestGroup1 = largestGroupSize1_1[i];
			}
		}
		else
		{
			largestGroupSize1_1[i] = 0;
		}
	}

    largestGroupSize0_1[nLoci - 1] = GET_BIT(nLoci, _c0, nLoci - 1) == GET_BIT(nLoci, target, nLoci - 1) ? 1 : 0;
    largestGroupSize1_0[nLoci - 1] = GET_BIT(nLoci, _c1, nLoci - 1) == GET_BIT(nLoci, target, nLoci - 1) ? 1 : 0;
	for (int i = nLoci - 2; i >= 0; i--)
	{
        if (GET_BIT(nLoci, _c1, i) == GET_BIT(nLoci, target, i))
		{
			largestGroupSize1_0[i] += largestGroupSize1_0[i + 1] + 1;
		}
		else
		{
			largestGroupSize1_0[i] = 0;
		}

        if (GET_BIT(nLoci, _c0, i) == GET_BIT(nLoci, target, i))
		{
			largestGroupSize0_1[i] += largestGroupSize0_1[i + 1] + 1;
		}
		else
		{
			largestGroupSize0_1[i] = 0;
		}
	}

	/* Determine best crossover point */
	std::vector<std::pair<int, bool> > bestCrossOverPoint;
	int largestGroup = std::max(largestGroup0, largestGroup1);

	for (int i = 0; i < nLoci - 1; i++)
	{
		if ((largestGroupSize0_0[i] + largestGroupSize1_0[i + 1]) > largestGroup)
		{
			largestGroup = largestGroupSize0_0[i] + largestGroupSize1_0[i + 1];
			bestCrossOverPoint.clear();
			bestCrossOverPoint.push_back(std::pair<int, bool>(i, true));
		}
		
		if ((largestGroupSize1_1[i] + largestGroupSize0_1[i + 1]) > largestGroup)
		{
			largestGroup = largestGroupSize1_1[i] + largestGroupSize0_1[i + 1];
			bestCrossOverPoint.clear();
			bestCrossOverPoint.push_back(std::pair<int, bool>(i, false));
		}
	}

	/* now generate the mask */
	if (bestCrossOverPoint.empty() && largestGroup0 > largestGroup1)
	{
		return std::vector<int>(_c0);
	}
	else if (bestCrossOverPoint.empty()) // largestGroup0 <= largestGroup1
	{
		return std::vector<int>(_c1);
	}
	else
	{
		std::vector<int> res;
		for (std::vector<std::pair<int, bool> >::const_iterator it = bestCrossOverPoint.begin(); 
			it != bestCrossOverPoint.end(); it++)
		{
			int val1 = 0;
			int val2 = 0;
			for (int i = 0; i < nLoci; i++)
			{
				if (i <= it->first)
				{
					val1 |= (_c0 & (1 << i));
					val2 |= (_c1 & (1 << i));
				}
				else
				{
					val1 |= (_c1 & (1 << i));
					val2 |= (_c0 & (1 << i));
				}
			}
			
			if (it->second)
			{
				res.push_back(val1);
			}
			else
			{
				res.push_back(val2);
			}
		}

		return res;
	}
}
