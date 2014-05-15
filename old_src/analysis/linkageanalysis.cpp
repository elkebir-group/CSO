/*
 * linkageanalysis.cpp
 *
 *  Created on: 15-apr-2009
 *      Author: M. El-Kebir
 */

#include "linkageanalysis.h"

LinkageAnalysis::LinkageAnalysis(const Data* pData)
  : _nLoci(pData->getNumberOfLoci())
{
  const Genotype& ideotype = pData->getIdeotype();
  const GenotypeSet& parents = pData->getParents();
  const DoubleMatrix& RM = pData->getRM();
  const double gamma = pData->getGamma();

  analyse(_nLoci, parents, RM, gamma, ideotype.getC0());
  update(_nLoci, ideotype.getC0());

  if (!ideotype.isHomozygous())
  {
    analyse(_nLoci, parents, RM, gamma, ideotype.getC1());
    update(_nLoci, ideotype.getC1());
  }
}

LinkageAnalysis::LinkageAnalysis(const int nLoci, const GenotypeSet& parents, 
  const Genotype& ideotype, const DoubleMatrix& RM, const double gamma)
  : _nLoci(nLoci)
{
  analyse(nLoci, parents, RM, gamma, ideotype.getC0());
  update(nLoci, ideotype.getC0());

  if (!ideotype.isHomozygous())
  {
    analyse(nLoci, parents, RM, gamma, ideotype.getC1());
    update(nLoci, ideotype.getC1());
  }
}

LinkageAnalysis::~LinkageAnalysis()
{
}

void LinkageAnalysis::analyse(const int nLoci, const GenotypeSet& parents, 
  const DoubleMatrix& RM, const double gamma, int targetChromosome)
{
  LinkageMatrix& matrix = _linkageMatrixMap[targetChromosome];
  matrix = LinkageMatrix(nLoci);
  for (int i = 0; i < nLoci; i++)
  {
    matrix[i] = std::vector<LinkagePopPair>(nLoci, LinkagePopPair(Unlinked, ULONG_MAX));
  }

  for (int i = 0; i < nLoci; i++)
  {
    for (int j = i + 1; j < nLoci; j++)
    {
            int val_i = GET_BIT(nLoci, targetChromosome, i);
            int val_j = GET_BIT(nLoci, targetChromosome, j);

      for (GenotypeSet::const_iterator it = parents.begin(); it != parents.end(); it++)
      {
        const Genotype& parent = *it;
        if (parent(nLoci, 0, i) == val_i)
        {
          if (parent(nLoci, 0, j) == val_j)
          {
            // no recombination
            unsigned long pop = probToPop(parent.computeProb(nLoci, RM, parent.getC0()), gamma);
            if (pop < matrix[i][j].second || matrix[i][j].first != Linked)
              matrix[i][j] = matrix[j][i] = LinkagePopPair(Linked, pop);
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
            unsigned long pop = probToPop(parent.computeProb(nLoci, RM, parent.getC1()), gamma);
            if (pop < matrix[i][j].second || matrix[i][j].first != Linked)
              matrix[i][j] = matrix[j][i] = LinkagePopPair(Linked, pop);
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

inline void LinkageAnalysis::update(const int nLoci, const int targetChromosome)
{
  // initialize lists
  _linkedMap[targetChromosome];
  _weaklyLinkedMap[targetChromosome];
  _unlinkedMap[targetChromosome];

  const LinkageMatrix& matrix = _linkageMatrixMap[targetChromosome];
  for (int i = 0; i < nLoci; i++)
  {
    for (int j = i + 1; j < nLoci; j++)
    {
      const LinkagePopPair& linkagePopPair = matrix[i][j];
      switch (linkagePopPair.first)
      {
      case Linked:
        _linkedMap[targetChromosome].push_back(LociPair(i, j));
        break;
      case WeaklyLinked:
        _weaklyLinkedMap[targetChromosome].push_back(LociPair(i, j));
        break;
      case Unlinked:
        _unlinkedMap[targetChromosome].push_back(LociPair(i, j));
        break;
      default:
        assert(false);
      }
    }
  }
}

LinkagePopPair LinkageAnalysis::getLinkagePopPair(int targetChromosome, int i, int j) const
{
  assert(0 <= i && i < Data::getInstance()->getNumberOfLoci());
  assert(0 <= j && j < Data::getInstance()->getNumberOfLoci());
  assert(_linkageMatrixMap.find(targetChromosome) != _linkageMatrixMap.end());

  return _linkageMatrixMap.find(targetChromosome)->second[i][j];
}

unsigned long LinkageAnalysis::getPopMax() const
{
  unsigned long max = 0;

  for (LinkageMatrixMap::const_iterator it = _linkageMatrixMap.begin();
    it != _linkageMatrixMap.end();
    it++)
  {
    for (int i = 0; i < _nLoci; i++)
    {
      for (int j = i + 1; j < _nLoci; j++)
      {
        if (max < it->second[i][j].second)
        {
          max = it->second[i][j].second;
        }
      }
    }
  }

  return max;
}
