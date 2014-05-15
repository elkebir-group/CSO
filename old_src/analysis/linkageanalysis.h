/*
 * linkageanalysis.h
 *
 *  Created on: 15-apr-2009
 *      Author: M. El-Kebir
 */

#ifndef LINKAGEANALYSIS_H_
#define LINKAGEANALYSIS_H_

#include "cso.h"
#include "data.h"
#include <vector>
#include <map>
#include "genotype.h"

typedef std::pair<LinkageType, unsigned long> LinkagePopPair;
typedef std::vector<std::vector<LinkagePopPair> > LinkageMatrix;

typedef std::pair<int, int> LociPair;
typedef std::list<LociPair> LociPairList;

typedef std::map<int, LociPairList> LociPairListMap;
typedef std::map<int, LinkageMatrix> LinkageMatrixMap;

/* PRE-CONDITION: ideotype can be obtained from parents */
class LinkageAnalysis
{
private:
  const int _nLoci;
  LinkageMatrixMap _linkageMatrixMap;
  LociPairListMap _linkedMap;
  LociPairListMap _weaklyLinkedMap;
  LociPairListMap _unlinkedMap;

  void analyse(const int nLoci, const GenotypeSet& parents,
    const DoubleMatrix& RM, const double gamma, int targetChromosome);
  void update(const int nLoci, const int targetChromosome);

public:
  LinkageAnalysis(const Data* pData);
  LinkageAnalysis(const int nLoci, const GenotypeSet& parents,
    const Genotype& ideotype, const DoubleMatrix& RM, const double gamma);
  ~LinkageAnalysis();
  LinkagePopPair getLinkagePopPair(int targetChromosome, int i, int j) const;
  LinkagePopPair getLinkagePopPair(int targetChromosome, const LociPair& lociPair) const;
  LinkageType getLinkageType(int targetChromosome, int i, int j) const;
  const LociPairList& getLinkedLoci(int targetChromosome) const;
  const LociPairList& getWeaklyLinkedLoci(int targetChromosome) const;
  const LociPairList& getUnlinkedLoci(int targetChromosome) const;
  unsigned long getPopMax() const;
};


inline const LociPairList& LinkageAnalysis::getLinkedLoci(int targetChromosome) const
{
  assert(_linkedMap.find(targetChromosome) != _linkedMap.end());
  return _linkedMap.find(targetChromosome)->second;
}

inline const LociPairList& LinkageAnalysis::getWeaklyLinkedLoci(int targetChromosome) const
{
  assert(_weaklyLinkedMap.find(targetChromosome) != _weaklyLinkedMap.end());
  return _weaklyLinkedMap.find(targetChromosome)->second;
}

inline const LociPairList& LinkageAnalysis::getUnlinkedLoci(int targetChromosome) const
{
  assert(_unlinkedMap.find(targetChromosome) != _unlinkedMap.end());
  return _unlinkedMap.find(targetChromosome)->second;
}

inline LinkagePopPair LinkageAnalysis::getLinkagePopPair(int targetChromosome, const LociPair& lociPair) const
{
  return getLinkagePopPair(targetChromosome, lociPair.first, lociPair.second);
}

inline LinkageType LinkageAnalysis::getLinkageType(int targetChromosome, int i, int j) const
{
  return getLinkagePopPair(targetChromosome, i, j).first;
}

#endif
