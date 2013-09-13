/*
 * dpitem.cpp
 *
 *  Created on: 19-feb-2009
 *      Author: s030858
 */

#include "dpitem.h"
#include "data.h"
#include <algorithm>

DpItem::DpItem()
  : GenotypeGamete(-1, -1)
  , _parent1(NULL)
  , _parent2(NULL)
  , _cumCross(INT_MAX)
  , _cumCrossover(INT_MAX)
  , _crossoverP1(INT_MAX)
  , _crossoverP2(INT_MAX)
  , _gen(INT_MAX)
  , _ancestors()
  , _updateFlag(false)
{
  // this constructor should not be called, but is needed for unordered_map
  assert(false);
  abort();
}

DpItem::DpItem(const Genotype& genotype, DpItem* parent1, DpItem* parent2, 
               double crossoverP1, double crossoverP2,
               unsigned long gen, unsigned long cumCrossover, 
               unsigned long cumCross, const GenotypeSet& ancestors, bool updateFlag)
  : GenotypeGamete(genotype)
  , _parent1(parent1)
  , _parent2(parent2)
  , _cumCross(cumCross)
  , _cumCrossover(cumCrossover)
  , _crossoverP1(crossoverP1)
  , _crossoverP2(crossoverP2)
  , _gen(gen)
  , _ancestors(ancestors)
  , _updateFlag(updateFlag)
{
}

DpItem::DpItem(const Genotype& genotype, bool updateFlag)
  : GenotypeGamete(genotype)
  , _parent1(NULL)
  , _parent2(NULL)
  , _cumCross(0)
  , _cumCrossover(0)
  , _crossoverP1(0)
  , _crossoverP2(0)
  , _gen(0)
  , _ancestors()
  , _updateFlag(updateFlag)
{
}

DpItem::~DpItem()
{
}

void DpItem::updateAttributesFast(const DpTable* pTable, bool updateFlag)
{
  if (_updateFlag == updateFlag)
  {
    // stop recursion, as this item has already been updated
    return;
  }

  _ancestors.clear();
  _cumCrossover = getCrossover();
  _cumCross = 1;

  if (!_parent1 && !_parent2)
  {
    // item is a leaf, stop the recursion
    _cumCross = 0;
    _gen = 0;

    // set updateFlag to new value
    _updateFlag = updateFlag;
    return;
  }

  //assert(_parent1 && _parent2);

  if (_parent1 && _parent1->_updateFlag != updateFlag)
  {
    // update attributes of parent1
    _parent1->updateAttributesFast(pTable, updateFlag);
  }
  if (_parent2 && _parent2->_updateFlag != updateFlag)
  {
    // update attributes of parent2
    _parent2->updateAttributesFast(pTable, updateFlag);
  }

  // new ancestor set is the union of the parental ancestor sets 
  // plus the two parents themselves
  if (_parent1 && _parent2)
  {
    set_union(_parent1->_ancestors.begin(), _parent1->_ancestors.end(),
      _parent2->_ancestors.begin(), _parent2->_ancestors.end(),
      inserter(_ancestors, _ancestors.begin()));
    _ancestors.insert(*_parent1);
    _ancestors.insert(*_parent2);
  }
  else if (_parent1)
  {
    _ancestors = _parent1->_ancestors;
    _ancestors.insert(*_parent1);
  }
  else
  {
    _ancestors = _parent2->_ancestors;
    _ancestors.insert(*_parent2);
  }


  // update _cumCrossover and _cumCross (number of internal nodes)
  for (GenotypeSet::const_iterator it = _ancestors.begin(); it != _ancestors.end(); it++)
  {
    const DpItem& item = pTable->getItem(*it);
    
    _cumCrossover += item.getCrossover();
    
    if (!item.isParent())
      _cumCross++;
  }

  // update _gen
  if (_parent1 && _parent2)
    _gen = 1 + std::max(_parent1->_gen, _parent2->_gen);
  else if (_parent1)
    _gen = 1 + _parent1->_gen;
  else
    _gen = 1 + _parent2->_gen;

  // set updateFlag to new value
  _updateFlag = updateFlag;
}

void DpItem::printItem(int nLoci, std::ostream& out, bool newLine) const
{
  printGenotype(nLoci, false);
  out << "\t";
  if (_parent1)
  {
    _parent1->printGenotype(nLoci, false);
  }
  else
  {
    for (int i = 0; i < nLoci; i++)
    {
      out << '-';
    }
    
    out << '/';
    
    for (int i = 0; i < nLoci; i++)
    {
      out << '-';
    }
  }
  out << '\t';
  if (_parent2)
  {
    _parent2->printGenotype(nLoci, false);
  }
  else
  {
    for (int i = 0; i < nLoci; i++)
    {
      out << '-';
    }
    
    out << '/';
    
    for (int i = 0; i < nLoci; i++)
    {
      out << '-';
    }
  }
  out << "\t" << _gen << "\t" << _cumCross << "\t" << _cumCrossover << "\t" 
    << Data::getInstance()->getCost(_cumCrossover, _gen, _cumCross);
  
  if (newLine) 
    out << std::endl;
}

void DpItem::printNode(int nLoci, const DpItem* pItem, std::ostream& out) const
{
  static const Genotype& ideotype = Data::getInstance()->getIdeotype();

  out << "\t\t";
  pItem->printGenotype(nLoci, false, out, "");
  
  out << " [fontcolor=white,style=rounded,shape=box,";
  out << "label=";

  out << "<<TABLE BORDER=\"0\" CELLSPACING=\"0\" CELLPADDING=\"4\" CELLBORDER=\"0\">" << std::endl;
  out << "\t\t\t<TR>" << std::endl;

  char buf[33];
  toBitstring(pItem->getC0(), nLoci, buf);
  for (int i = 0; i < nLoci; i++)
  {
    out << "\t\t\t\t<TD BGCOLOR=\"";

    if (ideotype(0, nLoci - i - 1) == (*pItem)(0, nLoci - i - 1))
      out << "blue";
    else 
      out << "red";

    out << "\">" << buf[i] << "</TD>" << std::endl;
  }

  out << "\t\t\t\t<TD WIDTH=\"40\" ROWSPAN=\"2\"><FONT COLOR=\"black\">";
  out << pItem->_cumCrossover << "<BR/>" 
    << pItem->_gen << "/" << pItem->_cumCross << "<BR/>" 
    << Data::getInstance()->getCost(pItem->_cumCrossover, pItem->_gen, pItem->_cumCross);
  out << "</FONT></TD>" << std::endl;

  out << "\t\t\t</TR>" << std::endl;

  out << "\t\t\t<TR>" << std::endl;

  toBitstring(pItem->getC1(), nLoci, buf);
  for (int i = 0; i < nLoci; i++)
  {
    out << "\t\t\t\t<TD BGCOLOR=\"";

    if (ideotype(1, nLoci - i - 1) == (*pItem)(1, nLoci - i - 1))
      out << "blue";
    else 
      out << "red";

    out << "\">" << buf[i] << "</TD>" << std::endl;
  }

  out << "\t\t\t</TR>" << std::endl;
  out << "\t\t</TABLE>>" << "];" << std::endl;
}

void DpItem::printNodes(const DpTable* pTable, int nLoci, std::ostream& out) const
{
  std::vector<DpItemPointerList> genotypesPerGeneration(_gen+1);
  
  for (GenotypeSet::const_iterator it = _ancestors.begin(); it != _ancestors.end(); it++)
  {
    const DpItem& item = pTable->getItem(*it);
    genotypesPerGeneration[item._gen].push_back(&item);
  }
  genotypesPerGeneration[_gen].push_back(this);

  for (unsigned long i = 0; i <= _gen; i++)
  {
    const DpItemPointerList& list = genotypesPerGeneration[i];
    
    if (list.size() != 0)
    {
      out << "\t{\n\t\trank = same;\n";
    }

    for (DpItemPointerList::const_iterator it = list.begin(); it != list.end(); it++)
    {
      printNode(nLoci, *it, out);

      //out << " [shape=Mrecord,";

      /*out << "label=\"{ {";

      char buf[33];
      toBitstring((*it)->getC0(), nLoci, buf);
      for (int i = 0; i < nLoci; i++)
      {
        out << buf[i];
        if (i < nLoci - 1) out << '|';
      }
      out << "} | {";
      
      toBitstring((*it)->getC1(), nLoci, buf);
      for (int i = 0; i < nLoci; i++)
      {
        out << buf[i];
        if (i < nLoci - 1) out << '|';
      }
      out << "} } | " << (*it)->_cumCrossover << "\\n" 
        << (*it)->_gen << "/" << (*it)->_cumCross << "\\n" 
        << Solver::_DP->getScore(*(*it)) << "\"];" << std::endl;
      */
    }

    if (list.size() != 0)
    {
      out << "\t}\n";
    }
  }
}

void DpItem::printEdges(std::ostream& out, int nLoci, GenotypeSet& printedGenotypes) const
{
  if (printedGenotypes.find(*this) == printedGenotypes.end())
  {
    //if (_parent1 != _parent2)
    //{
      printedGenotypes.insert(*this);

      if (_parent1)
      {
        out << '\t';
        _parent1->printGenotype(nLoci, false, out, "");
        out << " -> ";
        printGenotype(nLoci, false, out, "");
        out << " [label=\"" << _crossoverP1 << "\"]";
        out << ';' << std::endl;
        _parent1->printEdges(out, nLoci, printedGenotypes);
      }

      if (_parent2)
      {
        out << '\t';
        _parent2->printGenotype(nLoci, false, out, "");
        out << " -> ";
        printGenotype(nLoci, false, out, "");
        out << " [label=\"" << _crossoverP2 << "\"]";
        out << ';' << std::endl;
        _parent2->printEdges(out, nLoci, printedGenotypes);
      }
    }
    //if (_parent1 == _parent2 && _parent1)
    //{
    //  out << '\t';
    //  _parent1->printGenotype(nLoci, false, out, "");
    //  out << " -> ";
    //  printGenotype(nLoci, false, out, "");
    //  out << " [label=\"" << getCrossover() << "\"]";
    //  out << ';' << std::endl;

    //  printedGenotypes.insert(*this);
    //  _parent1->printEdges(out, nLoci, printedGenotypes);
    //}
  //}
}
