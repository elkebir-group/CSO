/*
 * dpitem.h
 *
 *  Created on: 19-feb-2009
 *      Author: s030858
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
  double _cumCrossover;
  double _crossoverP1;
  double _crossoverP2;
  unsigned long _gen;
  GenotypeSet _ancestors;
  bool _updateFlag;
  void printNode(int nLoci, const DpItem* pItem, std::ostream& out) const;

public:
  DpItem();
  DpItem(const Genotype& genotype, DpItem* parent1, DpItem* parent2, 
    double crossoverP1, double crossoverP2,
    unsigned long gen, unsigned long cumCrossover, 
    unsigned long cumCross, const GenotypeSet& ancestors, bool updateFlag);
  DpItem(const Genotype& genotype, bool updateFlag);
  ~DpItem();
  const DpItem* getParent1() const;
  const DpItem* getParent2() const;
  const GenotypeSet& getAncestors() const;
  bool isParent() const;
  void printItem(int nLoci, std::ostream& out = std::cout, bool newLine = true) const;
  void printEdges(std::ostream& out, int nLoci, GenotypeSet& printedGenotypes) const;
  void printNodes(const DpTable* pTable, int nLoci, std::ostream& out) const;
  double getCumCrossover() const;
  unsigned long getCumCross() const;
  unsigned long getGen() const;
  double getCrossover() const;
  double getCrossoverP1() const;
  double getCrossoverP2() const;
  void update(const DpItem& item);
  bool hasAncestor(const Genotype& genotype) const;
  void updateAttributesFast(const DpTable* pTable, bool updateFlag);
  bool getUpdateFlag() const;
};

inline double DpItem::getCrossoverP1() const
{
  return _crossoverP1;
}

inline double DpItem::getCrossoverP2() const
{
  return _crossoverP2;
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

inline double DpItem::getCumCrossover() const
{
  return _cumCrossover;
}

inline unsigned long DpItem::getGen() const
{
  return _gen;
}

inline double DpItem::getCrossover() const
{
  return _crossoverP1 + _crossoverP2;
}

inline void DpItem::update(const DpItem& item)
{
  assert(this != item._parent1 && this != item._parent2);
  _parent1 = item._parent1;
  _parent2 = item._parent2;
  _cumCross = item._cumCross;
  _cumCrossover = item._cumCrossover;
  _crossoverP1 = item._crossoverP1;
  _crossoverP2 = item._crossoverP2;
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
