/*
 * dpitem.h
 *
 *  Created on: 19-feb-2009
 *      Author: M. El-Kebir
 */

#ifndef DPITEM_H_
#define DPITEM_H_

#include "csodp.h"
#include "genotypegamete.h"
#include "dptable.h"

class DpItem : public GenotypeGamete
{
private:
	DpItem* _parent1;
	DpItem* _parent2;
	unsigned long _cumCross;
	unsigned long _cumPop;
	unsigned long _pop;
	unsigned long _gen;
	GenotypeSet _ancestors;
	bool _updateFlag;
	void printNode(int nLoci, const DpItem* pItem, std::ostream& out) const;

public:
	DpItem();
	DpItem(const Genotype& genotype, DpItem* parent1, DpItem* parent2, 
		unsigned long pop, unsigned long gen, unsigned long cumPop, 
		unsigned long cumCross, const GenotypeSet& ancestors, bool updateFlag);
	DpItem(const Genotype& genotype, DpItem* parent1, DpItem* parent2, 
		unsigned long pop, unsigned long gen, const GenotypeSet& ancestors, bool updateFlag);
	DpItem(const Genotype& genotype, bool updateFlag);
	~DpItem();
	const DpItem* getParent1() const;
	const DpItem* getParent2() const;
	const GenotypeSet& getAncestors() const;
	bool isParent() const;
	void printItem(int nLoci, std::ostream& out = std::cout, bool newLine = true) const;
	void printEdges(std::ostream& out, int nLoci, GenotypeSet& printedGenotypes) const;
	void printNodes(const DpTable* pTable, int nLoci, std::ostream& out) const;
	unsigned long getCumPop() const;
	unsigned long getCumCross() const;
	unsigned long getGen() const;
	unsigned long getPop() const;
	void update(const DpItem& item);
	bool hasAncestor(const Genotype& genotype) const;
	void updateAttributes(const DpTable* pTable, bool updateFlag, DpItemSet& newItemSet);
	void updateAttributesFast(const DpTable* pTable, bool updateFlag);
	bool getUpdateFlag() const;
};

inline DpItem::DpItem()
	: GenotypeGamete(-1, -1)
	, _parent1(NULL)
	, _parent2(NULL)
	, _cumCross(INT_MAX)
	, _cumPop(INT_MAX)
	, _pop(INT_MAX)
	, _gen(INT_MAX)
	, _ancestors()
	, _updateFlag(false)
{
	// this constructor should not be called, but is needed for unordered_map
	assert(false);
	abort();
}

inline DpItem::DpItem(const Genotype& genotype, DpItem* parent1, DpItem* parent2, 
	unsigned long pop, unsigned long gen, unsigned long cumPop, 
	unsigned long cumCross, const GenotypeSet& ancestors, bool updateFlag)
	: GenotypeGamete(genotype)
	, _parent1(parent1)
	, _parent2(parent2)
	, _cumCross(cumCross)
	, _cumPop(cumPop)
	, _pop(pop)
	, _gen(gen)
	, _ancestors(ancestors)
	, _updateFlag(updateFlag)
{
}

inline DpItem::DpItem(const Genotype& genotype, bool updateFlag)
	: GenotypeGamete(genotype)
	, _parent1(NULL)
	, _parent2(NULL)
	, _cumCross(0)
	, _cumPop(1)
	, _pop(1)
	, _gen(0)
	, _ancestors()
	, _updateFlag(updateFlag)
{
}

inline DpItem::~DpItem()
{
}

inline const DpItem* DpItem::getParent1() const
{
	return _parent1;
}

inline const DpItem* DpItem::getParent2() const
{
	return _parent2;
}

inline bool DpItem::isParent() const
{
	return _parent1 == NULL || _parent2 == NULL;
}

inline const GenotypeSet& DpItem::getAncestors() const
{
	return _ancestors;
}

inline unsigned long DpItem::getCumCross() const
{
	return _cumCross;
}

inline unsigned long DpItem::getCumPop() const
{
	return _cumPop;
}

inline unsigned long DpItem::getGen() const
{
	return _gen;
}

inline unsigned long DpItem::getPop() const
{
	return _pop;
}

inline void DpItem::update(const DpItem& item)
{
	assert(this != item._parent1 && this != item._parent2);
	_parent1 = item._parent1;
	_parent2 = item._parent2;
	_cumCross = item._cumCross;
	_cumPop = item._cumPop;
	_pop = item._pop;
	_gen = item._gen;
	_ancestors = item._ancestors;
	_updateFlag = item._updateFlag;
}

inline bool DpItem::hasAncestor(const Genotype& genotype) const
{
	return _ancestors.find(genotype) != _ancestors.end();
}

inline bool DpItem::getUpdateFlag() const
{
	return _updateFlag;
}

#endif /* DPITEM_H_ */
