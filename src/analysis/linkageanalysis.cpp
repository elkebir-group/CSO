/*
 * linkageanalysis.cpp
 *
 *  Created on: 15-apr-2009
 *      Author: s030858
 */

#include "linkageanalysis.h"

LinkageAnalysis::LinkageAnalysis(const Data* pData)
	: _linkageMatrixC0(pData->getNumberOfLoci())
	, _linkageMatrixC1(pData->getNumberOfLoci())
{
	const int nLoci = pData->getNumberOfLoci();
	const Genotype& ideotype = pData->getIdeotype();
	const GenotypeSet& parents = pData->getParents();
	const DoubleMatrix& RM = pData->getRM();
	const double gamma = pData->getGamma();

	for (int i = 0; i < nLoci; i++)
	{
		_linkageMatrixC0[i] = _linkageMatrixC1[i] = 
			std::vector<LinkagePopPair>(nLoci, LinkagePopPair(Unlinked, ULONG_MAX));
	}

	analyse(nLoci, parents, RM, gamma, _linkageMatrixC0, ideotype.getC0());
	analyse(nLoci, parents, RM, gamma, _linkageMatrixC1, ideotype.getC1());

	for (int i = 0; i < nLoci; i++)
	{
		for (int j = i + 1; j < nLoci; j++)
		{
			update(_linkageMatrixC0, LociPair(i, j), _linkedC0, _weaklyLinkedC0, _unlinkedC0);
			update(_linkageMatrixC1, LociPair(i, j), _linkedC1, _weaklyLinkedC1, _unlinkedC1);
		}
	}
}

LinkageAnalysis::LinkageAnalysis(const int nLoci, const GenotypeSet& parents, 
	const Genotype& ideotype, const DoubleMatrix& RM, const double gamma)
	: _linkageMatrixC0(nLoci)
	, _linkageMatrixC1(nLoci)
{
	for (int i = 0; i < nLoci; i++)
	{
		_linkageMatrixC0[i] = _linkageMatrixC1[i] = 
			std::vector<LinkagePopPair>(nLoci, LinkagePopPair(Unlinked, ULONG_MAX));
	}

	analyse(nLoci, parents, RM, gamma, _linkageMatrixC0, ideotype.getC0());
	analyse(nLoci, parents, RM, gamma, _linkageMatrixC1, ideotype.getC1());

	for (int i = 0; i < nLoci; i++)
	{
		for (int j = i + 1; j < nLoci; j++)
		{
			update(_linkageMatrixC0, LociPair(i, j), _linkedC0, _weaklyLinkedC0, _unlinkedC0);
			update(_linkageMatrixC1, LociPair(i, j), _linkedC1, _weaklyLinkedC1, _unlinkedC1);
		}
	}
}

LinkageAnalysis::~LinkageAnalysis()
{
}

void LinkageAnalysis::analyse(const int nLoci, const GenotypeSet& parents, 
	const DoubleMatrix& RM, const double gamma, LinkageMatrix& matrix, int target)
{
	for (int i = 0; i < nLoci; i++)
	{
		for (int j = i + 1; j < nLoci; j++)
		{
      int val_i = GET_BIT(nLoci, target, i);
      int val_j = GET_BIT(nLoci, target, j);

			for (GenotypeSet::const_iterator it = parents.begin(); it != parents.end(); it++)
			{
				const Genotype& parent = *it;
        if (parent(nLoci, 0, i) == val_i)
				{
          if (parent(nLoci, 0, j) == val_j)
					{
						// no recombination
						matrix[i][j] = matrix[j][i] = LinkagePopPair(Linked, probToPop(0.5 * (1 - RM[i][j]), gamma));
					}
					else if (parent(nLoci, 1, j) == val_j)
					{
						// recombination needed
						unsigned long pop = probToPop(0.5 * RM[i][j], gamma);
						if (pop < matrix[i][j].second || matrix[i][j].first == Unlinked)
							matrix[i][j] = matrix[j][i] = LinkagePopPair(WeaklyLinked, pop);
					}
					else
					{
						// recombination needed
						unsigned long pop = probToPop(0.5 * RM[i][j], gamma);
						if (pop < matrix[i][j].second)
							matrix[i][j] = matrix[j][i] = LinkagePopPair(Unlinked, pop);
					}
				}
				if (parent(nLoci, 1, i) == val_i)
				{
					if (parent(nLoci, 1, j) == val_j)
					{
						// no recombination
						matrix[i][j] = matrix[j][i] = LinkagePopPair(Linked, probToPop(0.5 * (1 - RM[i][j]), gamma));
					}
					else if (parent(nLoci, 0, j) == val_j)
					{
						// recombination needed
						unsigned long pop = probToPop(0.5 * RM[i][j], gamma);
						if (pop < matrix[i][j].second || matrix[i][j].first == Unlinked)
							matrix[i][j] = matrix[j][i] = LinkagePopPair(WeaklyLinked, pop);
					}
					else
					{
						// recombination needed
						unsigned long pop = probToPop(0.5 * RM[i][j], gamma);
						if (pop < matrix[i][j].second)
							matrix[i][j] = matrix[j][i] = LinkagePopPair(Unlinked, pop);
					}
				}
			}
		}
	}
}

LinkagePopPair LinkageAnalysis::getLinkagePopPair(bool useC0, int i, int j) const
{
	assert(0 <= i && i < Data::getInstance()->getNumberOfLoci());
	assert(0 <= j && j < Data::getInstance()->getNumberOfLoci());

	if (useC0)
	{
		return _linkageMatrixC0[i][j];
	}
	else
	{
		return _linkageMatrixC1[i][j];
	}
}

unsigned long LinkageAnalysis::getPopMax() const
{
	int nLoci = (int) _linkageMatrixC0.size();
	unsigned long res = 0;

	for (int i = 0; i < nLoci; i++)
	{
		for (int j = i + 1; j < nLoci; j++)
		{
			if (res < _linkageMatrixC0[i][j].second) 
				res = _linkageMatrixC0[i][j].second;
			
			if (res < _linkageMatrixC1[i][j].second) 
				res = _linkageMatrixC1[i][j].second;
		}
	}

	return res;
}
