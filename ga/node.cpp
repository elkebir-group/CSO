/*
 * node.cpp
 *
 *  Created on: 24-mar-2009
 *      Author: M. El-Kebir
 */

#include "node.h"
#include "leafnode.h"
#include "innernode.h"
#include "genotypetable.h"
#include "data.h"

NodePtr Node::clone(const Node* pNode)
{
	const InnerNode* pInnerNode = dynamic_cast<const InnerNode*>(pNode);
	const LeafNode* pLeafNode = dynamic_cast<const LeafNode*>(pNode);

	if (pInnerNode)
	{
		return InnerNode::clone(pInnerNode);
	}
	else
	{
		assert(pLeafNode);
		return LeafNode::clone(pLeafNode);
	}
}

NodePtr Node::clone(const NodePtr pNode)
{
	return Node::clone(pNode.get());
}

NodePtr Node::createSplitCross(const GenotypeSet& parents, int targetGamete)
{
	static int nLoci = Data::getInstance()->getNumberOfLoci();
	int numberOfParents = (int) parents.size();
	if (numberOfParents == 1)
	{
		const GenotypeGamete& parent = GenotypeTable::getInstance()->getGenotype(*parents.begin());
		//return LeafNode::create(parent, 1.0); // performt veel beter ;)
		return LeafNode::create(parent, randDouble() < 0.5, parent.getGameteMostAlike(nLoci, targetGamete)._c, 1.0);
	}

	GenotypeSet::const_iterator it = parents.begin();
	GenotypeSet parents1, parents2;
	int splitPoint = randInt(numberOfParents - 1);
	for (int i = 0; i <= splitPoint; i++)
	{
		parents1.insert(*it);
		it++;
	}
	for (int i = splitPoint; i < numberOfParents - 1; i++)
	{
		parents2.insert(*it);
		it++;
	}

	NodePtr p1 = createSplitCross(parents1, targetGamete);
	NodePtr p2 = createSplitCross(parents2, targetGamete);

	int p1Gamete = p1->getGamete();
	int p2Gamete = p2->getGamete();

	const GenotypeGamete& genotypeChild = GenotypeTable::getInstance()->getGenotype(p1Gamete, p2Gamete);
	int gameteChild = genotypeChild.getGameteMostAlike(nLoci, targetGamete)._c;

	if (gameteChild == p1Gamete)
	{
		return p1;
	}
	else if (gameteChild == p2Gamete)
	{
		return p2;
	}

	return InnerNode::create(p1, p2, randInt(2) == 1, (BackcrossType) randInt(3), gameteChild);
}

void Node::printNode(std::ostream& out, const Genotype& genotype, unsigned long gen, unsigned long pop, unsigned long cross)
{
	static const Genotype& ideotype = Data::getInstance()->getIdeotype();
	static const int nLoci = Data::getInstance()->getNumberOfLoci();

	out << "\t\t";
	genotype.printGenotype(nLoci, false, out, "");
	
	out << " [fontcolor=white,style=rounded,shape=box,";
	out << "label=";

	out << "<<TABLE BORDER=\"0\" CELLSPACING=\"0\" CELLPADDING=\"4\" CELLBORDER=\"0\">" << std::endl;
	out << "\t\t\t<TR>" << std::endl;

	char buf[33];
	toBitstring(genotype.getC0(), nLoci, buf);
	for (int i = 0; i < nLoci; i++)
	{
		out << "\t\t\t\t<TD BGCOLOR=\"";

		if (ideotype(0, nLoci - i - 1) == genotype(0, nLoci - i - 1))
			out << "blue";
		else 
			out << "red";

		out << "\">" << buf[i] << "</TD>" << std::endl;
	}

	out << "\t\t\t\t<TD WIDTH=\"40\" ROWSPAN=\"2\"><FONT COLOR=\"black\">";
	out << pop << "<BR/>" 
		<< gen << "/" << cross << "<BR/>" 
		<< Data::getInstance()->getCost(pop, gen, cross);
	out << "</FONT></TD>" << std::endl;

	out << "\t\t\t</TR>" << std::endl;

	out << "\t\t\t<TR>" << std::endl;

	toBitstring(genotype.getC1(), nLoci, buf);
	for (int i = 0; i < nLoci; i++)
	{
		out << "\t\t\t\t<TD BGCOLOR=\"";

		if (ideotype(1, nLoci - i - 1) == genotype(1, nLoci - i - 1))
			out << "blue";
		else 
			out << "red";

		out << "\">" << buf[i] << "</TD>" << std::endl;
	}

	out << "\t\t\t</TR>" << std::endl;
	out << "\t\t</TABLE>>" << "];" << std::endl;
}

void Node::printEdge(std::ostream& out, const Genotype& parent1, const Genotype& parent2, const Genotype& child)
{
	int nLoci = Data::getInstance()->getNumberOfLoci();
	double gamma = Data::getInstance()->getGamma();
	const DoubleMatrix& RM = Data::getInstance()->getRM();

	unsigned int pop = child.computePop(nLoci, RM, gamma, parent1, parent2);

	out << '\t';
	parent1.printGenotype(nLoci, false, out, "");
	out << " -> ";
	child.printGenotype(nLoci, false, out, "");
	out << " [label=\"" << pop << "\"]";
	out << ';' << std::endl;

	if (parent1 != parent2)
	{
		out << '\t';
		parent2.printGenotype(nLoci, false, out, "");
		out << " -> ";
		child.printGenotype(nLoci, false, out, "");
		out << " [label=\"" << pop << "\"]";
		out << ';' << std::endl;
	}
}
