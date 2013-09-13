/*
 * linkageanalysis.h
 *
 *  Created on: 15-apr-2009
 *      Author: s030858
 */

#ifndef LINKAGEANALYSIS_H_
#define LINKAGEANALYSIS_H_

#include "cso.h"
#include "data.h"
#include <vector>

typedef std::pair<LinkageType, unsigned long> LinkagePopPair;
typedef std::vector<std::vector<LinkagePopPair> > LinkageMatrix;
typedef std::pair<int, int> LociPair;
typedef std::list<LociPair> LociPairList;

/* PRE-CONDITION: ideotype can be obtained from parents */
class LinkageAnalysis
{
private:
	LinkageMatrix _linkageMatrixC0;
	LinkageMatrix _linkageMatrixC1;
	LociPairList _linkedC0;
	LociPairList _linkedC1;
	LociPairList _weaklyLinkedC0;
	LociPairList _weaklyLinkedC1;
	LociPairList _unlinkedC0;
	LociPairList _unlinkedC1;

	void analyse(const int nLoci, const GenotypeSet& parents, 
		const DoubleMatrix& RM, const double gamma, LinkageMatrix& matrix, int target);
	void update(LinkageMatrix& matrix, const LociPair& lociPair, 
		LociPairList& linked, LociPairList& weaklyLinked, LociPairList& unlinked);

public:
	LinkageAnalysis(const Data* pData);
	LinkageAnalysis(const int nLoci, const GenotypeSet& parents, 
		const Genotype& ideotype, const DoubleMatrix& RM, const double gamma);
	~LinkageAnalysis();
	LinkagePopPair getLinkagePopPair(bool useC0, int i, int j) const;
	LinkagePopPair getLinkagePopPair(bool useC0, const LociPair& lociPair) const;
	const LociPairList& getLinkedLoci(bool useC0) const;
	const LociPairList& getWeaklyLinkedLoci(bool useC0) const;
	const LociPairList& getUnlinkedLoci(bool useC0) const;
	unsigned long getPopMax() const;
};

inline const LociPairList& LinkageAnalysis::getLinkedLoci(bool useC0) const
{
	if (useC0)
		return _linkedC0;
	else
		return _linkedC1;
}

inline const LociPairList& LinkageAnalysis::getWeaklyLinkedLoci(bool useC0) const
{
	if (useC0)
		return _weaklyLinkedC0;
	else
		return _weaklyLinkedC1;
}

inline const LociPairList& LinkageAnalysis::getUnlinkedLoci(bool useC0) const
{
	if (useC0)
		return _unlinkedC0;
	else
		return _unlinkedC1;
}

inline LinkagePopPair LinkageAnalysis::getLinkagePopPair(bool useC0, const LociPair& lociPair) const
{
	return getLinkagePopPair(useC0, lociPair.first, lociPair.second);
}

inline void LinkageAnalysis::update(LinkageMatrix& matrix, const LociPair& lociPair, 
	LociPairList& linkedList, LociPairList& weaklyLinkedList, LociPairList& unlinkedList)
{
	const LinkagePopPair& linkagePopPair = matrix[lociPair.first][lociPair.second];
	switch (linkagePopPair.first)
	{
	case Linked:
		linkedList.push_back(lociPair);
		break;
	case WeaklyLinked:
		weaklyLinkedList.push_back(lociPair);
		break;
	case Unlinked:
		unlinkedList.push_back(lociPair);
		break;
	default:
		assert(false);
	}
}

#endif
