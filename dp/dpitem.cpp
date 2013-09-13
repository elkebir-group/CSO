/*
 * dpitem.cpp
 *
 *  Created on: 19-feb-2009
 *      Author: M. El-Kebir
 */

#include "solver.h"
#include "dpitem.h"
#include <algorithm>

void DpItem::updateAttributesFast(const DpTable* pTable, bool updateFlag)
{
	if (_updateFlag == updateFlag)
	{
		// stop recursion, as this item has already been updated
		return;
	}

	_ancestors.clear();
	_cumPop = _pop;
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

	assert(_parent1 && _parent2);

	if (_parent1->_updateFlag != updateFlag)
	{
		// update attributes of parent1
		_parent1->updateAttributesFast(pTable, updateFlag);
	}
	if (_parent2->_updateFlag != updateFlag)
	{
		// update attributes of parent2
		_parent2->updateAttributesFast(pTable, updateFlag);
	}

	// new ancestor set is the union of the parental ancestor sets 
	// plus the two parents themselves
	set_union(_parent1->_ancestors.begin(), _parent1->_ancestors.end(),
		_parent2->_ancestors.begin(), _parent2->_ancestors.end(),
		inserter(_ancestors, _ancestors.begin()));

	_ancestors.insert(*_parent1);
	_ancestors.insert(*_parent2);

	// update _cumPop and _cumCross (number of internal nodes)
	for (GenotypeSet::const_iterator it = _ancestors.begin(); it != _ancestors.end(); it++)
	{
		const DpItem& item = pTable->getItem(*it);
		
		_cumPop += item._pop;
		
		if (!item.isParent())
			_cumCross++;
	}

	// update _gen
	_gen = 1 + std::max(_parent1->_gen, _parent2->_gen);

	// set updateFlag to new value
	_updateFlag = updateFlag;
}

void DpItem::updateAttributes(const DpTable* pTable, bool updateFlag, DpItemSet& newItemSet)
{
	if (_updateFlag == updateFlag)
	{
		// stop recursion, as this item has already been updated
		return;
	}

	GenotypeSet newAncestors;
	unsigned long newCumPop = _pop;
	unsigned long newCumCross = 1;

	if (!_parent1 && !_parent2)
	{
		// item is a leaf, stop the recursion
		_gen = 0;
		_cumCross = 0;

		// set updateFlag to new value
		_updateFlag = updateFlag;
		return;
	}

	assert(_parent1 && _parent2);

	if (_parent1->_updateFlag != updateFlag)
	{
		// update attributes of parent1
		_parent1->updateAttributes(pTable, updateFlag, newItemSet);
	}
	if (_parent2->_updateFlag != updateFlag)
	{
		// update attributes of parent2
		_parent2->updateAttributes(pTable, updateFlag, newItemSet);
	}

	// new ancestor set is the union of the parental ancestor sets 
	// plus the two parents themselves
	set_union(_parent1->_ancestors.begin(), _parent1->_ancestors.end(),
		_parent2->_ancestors.begin(), _parent2->_ancestors.end(),
		inserter(newAncestors, newAncestors.begin()));

	newAncestors.insert(*_parent1);
	newAncestors.insert(*_parent2);

	// update _cumPop and _cumCross (number of internal nodes)
	for (GenotypeSet::const_iterator it = newAncestors.begin(); it != newAncestors.end(); it++)
	{
		const DpItem& item = pTable->getItem(*it);
		
		newCumPop += item.getPop();
		
		if (!item.isParent())
			newCumCross++;
	}

	// update _gen
	unsigned long newGen = 1 + std::max(_parent1->_gen, _parent2->_gen);

	if (_gen != newGen || _cumPop != newCumPop || _cumCross != newCumCross || newAncestors != _ancestors)
	{
		newItemSet.insert(*this);
	}
	
	_gen = newGen;
	_cumCross = newCumCross;
	_cumPop = newCumPop;
	_ancestors = newAncestors;

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
	out << "\t" << _gen << "\t" << _cumCross << "\t" << _cumPop << "\t" << Solver::_DP->getScore(*this);
	
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
	out << pItem->_cumPop << "<BR/>" 
		<< pItem->_gen << "/" << pItem->_cumCross << "<BR/>" 
		<< Solver::_DP->getScore(*pItem);
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
			out << "} } | " << (*it)->_cumPop << "\\n" 
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
		if (_parent1 && _parent2)
		{
			if (_parent1 != _parent2)
			{
				out << '\t';
				_parent1->printGenotype(nLoci, false, out, "");
				out << " -> ";
				printGenotype(nLoci, false, out, "");
				out << " [label=\"" << _pop << "\"]";
				out << ';' << std::endl;

				out << '\t';
				_parent2->printGenotype(nLoci, false, out, "");
				out << " -> ";
				printGenotype(nLoci, false, out, "");
				out << " [label=\"" << _pop << "\"]";
				out << ';' << std::endl;

				printedGenotypes.insert(*this);

				_parent1->printEdges(out, nLoci, printedGenotypes);
				_parent2->printEdges(out, nLoci, printedGenotypes);
			}
			else
			{
				out << '\t';
				_parent1->printGenotype(nLoci, false, out, "");
				out << " -> ";
				printGenotype(nLoci, false, out, "");
				out << " [label=\"" << _pop << "\"]";
				out << ';' << std::endl;

				printedGenotypes.insert(*this);
				_parent1->printEdges(out, nLoci, printedGenotypes);
			}
		}
	}
}
