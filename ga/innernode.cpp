/*
 * innernode.cpp
 *
 *  Created on: 20-mar-2009
 *      Author: M. El-Kebir
 */

#include "innernode.h"
#include "leafnode.h"
#include "genotypetable.h"
#include "data.h"

InnerNode::InnerNode(const InnerNode& innerNode)
	: Node(innerNode)
	, _genotype(innerNode._genotype)
	, _self(innerNode._self)
	, _backcross(innerNode._backcross)
	, _localPop(innerNode._localPop)
	, _gen(innerNode._gen)
{
	_p1 = Node::clone(innerNode._p1);
	_p2 = Node::clone(innerNode._p2);
}

InnerNode& InnerNode::operator =(const InnerNode& innerNode)
{
	if (this == &innerNode)
	{
		// self assignment
		return *this;
	}

	_genotype = innerNode._genotype;
	_self = innerNode._self;
	_backcross = innerNode._backcross;
	_localPop = innerNode._localPop;
	_gen = innerNode._gen;
	_gamete = innerNode._gamete;
	_p1 = Node::clone(innerNode._p1);
	_p2 = Node::clone(innerNode._p2);

	return *this;
}

InnerNode::InnerNode(NodePtr p1, NodePtr p2)
	: Node()
	, _genotype(p1->getGamete(), p2->getGamete())
	, _p1(p1)
	, _p2(p2)
	, _self(randInt(2) == 1)
	, _backcross((BackcrossType) randInt(3))
	, _localPop(0)
	, _gen(0)
{
	update(true);
}

InnerNode::InnerNode(NodePtr p1, NodePtr p2, bool self, BackcrossType backcross, int gamete)
	: Node(gamete)
	, _genotype(p1->getGamete(), p2->getGamete())
	, _p1(p1)
	, _p2(p2)
	, _self(self)
	, _backcross(backcross)
	, _localPop(0)
	, _gen(0)
{
	update(false);
}

InnerNode::~InnerNode()
{
}

void InnerNode::update(bool pickGamete)
{
	int nLoci = Data::getInstance()->getNumberOfLoci();
	double gamma = Data::getInstance()->getGamma();
	const DoubleMatrix& RM = Data::getInstance()->getRM();
	//int popMax = Data::getInstance()->getPopMax();

	const Genotype& genotypeP1 = _p1->getGenotype();
	const Genotype& genotypeP2 = _p2->getGenotype();
	
	int gameteP1 = _p1->getGamete();
	int gameteP2 = _p2->getGamete();

	if (pickGamete)
	{
		const GenotypeGamete& genotypeGamete = GenotypeTable::getInstance()->getGenotype(gameteP1, gameteP2);
		_gamete = genotypeGamete.getGameteProportionate()._c;
	}
	_localPop = _genotype.computePop(nLoci, RM, gamma, genotypeP1, genotypeP2);

	//if (_pGenotype.isHomozygous())
	//{
	//	_self = false;
	//}

	if (_backcross == BackcrossP1)
	{
		const Genotype newGenotype = Genotype(gameteP1, _gamete);
		unsigned long pop = newGenotype.computePop(nLoci, RM, gamma, _genotype, genotypeP1);
		//if (pop > popMax)
		//{
			// try backcross with p2
			//pNewGenotype = 
			//	GenotypeTable::getInstance()->getGenotype(gameteP2._c, _gamete._c);
			//
			//pop = pNewGenotype->computePop(nLoci, RM, gamma, *_pGenotype, genotypeP2);
			//if (pop > popMax)
			//{
			//	// backcross with p2 was also too expensive, no backcross
			//	_backcross = NoBackcross;
			//}
			//else
			//{
			//	_backcross = BackcrossP2;
			//	_localPop += pop;
			//	_pGenotype = pNewGenotype;
			//}
		//}
		//else
		{
			_localPop += pop;
			_genotype = newGenotype;
		}
	}
	else if (_backcross == BackcrossP2)
	{
		const Genotype newGenotype = Genotype(gameteP2, _gamete);
		
		unsigned long pop = newGenotype.computePop(nLoci, RM, gamma, _genotype, genotypeP2);
		//if (pop > popMax)
		//{
			// try backcross with p1
			//pNewGenotype = 
			//	GenotypeTable::getInstance()->getGenotype(gameteP1._c, _gamete._c);
			//
			//pop = pNewGenotype->computePop(nLoci, RM, gamma, *_pGenotype, genotypeP1);
			//if (pop > popMax)
			//{
				// backcross with p1 was also too expensive, no backcross
			//	_backcross = NoBackcross;
			//}
			//else
			//{
			//	_backcross = BackcrossP1;
			//	_localPop += pop;
			//	_pGenotype = pNewGenotype;
			//}
		//}
		//else
		{
			_localPop += pop;
			_genotype = newGenotype;
		}
	}

	if (_self)
	{
		//if (_pGenotype->isHomozygous())
		//{
		//	_self = false;
		//}
		//else
		{
			const Genotype newGenotype = Genotype(_gamete, _gamete);
			
			unsigned long pop = newGenotype.computePop(nLoci, RM, gamma, _genotype, _genotype);
			//if (pop > popMax)
			//{
			//	_self = false;
			//}
			//else
			{
				_localPop += pop;
				_genotype = newGenotype;
			}
		}
	}
}

double InnerNode::getRecc() const
{
	double reccP1 = _p1->getRecc();
	double reccP2 = _p2->getRecc();

	double reccRes = 0.5 * (reccP1 + reccP2);
	
	if (_backcross == BackcrossP1)
	{
		reccRes = 0.5 * (reccRes + reccP1);
	}
	else if (_backcross == BackcrossP2)
	{
		reccRes = 0.5 * (reccRes + reccP2);
	}

	return reccRes;
}

void InnerNode::getAllNodes(NodePtrSet& nodeSet) const
{
	if (nodeSet.find(NodePtr(_this)) == nodeSet.end())
	{
		nodeSet.insert(NodePtr(_this));
		_p1->getAllNodes(nodeSet);
		_p2->getAllNodes(nodeSet);
	}
}

NodePtr InnerNode::swapNodes(const NodePtr oldNode, const NodePtr newNode)
{
	if (NodePtr(_this) == oldNode)
	{
		return newNode;	
	}
	else
	{
		NodePtr newP1 = _p1->swapNodes(oldNode, newNode);
		NodePtr newP2 = _p1 == _p2 ? newP1 : _p2->swapNodes(oldNode, newNode);

		if (newP1->getGamete() == _p1->getGamete() &&
			newP2->getGamete() == _p2->getGamete())
		{
			return create(newP1, newP2, _self, _backcross, _gamete);
		}
		else
		{
			return create(newP1, newP2);
		}
	}
}

NodePtr InnerNode::create(NodePtr p1, NodePtr p2)
{
	InnerNode* pInnerNode = new InnerNode(p1, p2);
	NodePtr res(pInnerNode);
	pInnerNode->_this = res;
	return res;
}

NodePtr InnerNode::create(NodePtr p1, NodePtr p2, bool self, BackcrossType backcross, int gamete)
{
	InnerNode* pInnerNode = new InnerNode(p1, p2, self, backcross, gamete);
	NodePtr res(pInnerNode);
	pInnerNode->_this = res;
	return res;
}

NodePtr InnerNode::clone(const InnerNode* pInnerNode)
{
	InnerNode* pClone = new InnerNode(*pInnerNode);
	NodePtr res(pClone);
	pClone->_this = res;
	return res;
}

void InnerNode::updateGen()
{
	_p1->updateGen();
	if (_p1 != _p2)
		_p2->updateGen();
	
	_gen = getLocalGen() + std::max(_p1->getGen(), _p2->getGen());
}

void InnerNode::printEdges(std::ostream& out) const
{
	Genotype child1(_p1->getGamete(), _p2->getGamete());
	printEdge(out, _p1->getGenotype(), _p2->getGenotype(), child1);
	
	if (_backcross == BackcrossP1)
	{
		Genotype child2(_p1->getGamete(), _gamete);
		printEdge(out, _p1->getGenotype(), child1, child2);

		if (_self)
		{
			Genotype child3(_gamete, _gamete);
			printEdge(out, child2, child2, child3);
		}
	}
	else if (_backcross == BackcrossP2)
	{
		Genotype child2(_p2->getGamete(), _gamete);
		printEdge(out, _p2->getGenotype(), child1, child2);

		if (_self)
		{
			Genotype child3(_gamete, _gamete);
			printEdge(out, child2, child2, child3);
		}
	}
	else if (_self)
	{
		Genotype child2(_gamete, _gamete);
		printEdge(out, child1, child1, child2);
	}
}

void InnerNode::printNode(std::ostream& out, int level) const
{
	assert(0 <= level && level < (int) getLocalGen());

	int nLoci = Data::getInstance()->getNumberOfLoci();
	double gamma = Data::getInstance()->getGamma();
	const DoubleMatrix& RM = Data::getInstance()->getRM();

	unsigned long gen = _gen - getLocalGen();
	unsigned long pop = 0, cross = 0;

	NodePtrSet ancestors;
	getAllNodes(ancestors);
	for (NodePtrSet::const_iterator it = ancestors.begin();
		it != ancestors.end(); it++)
	{
		cross += (*it)->getLocalCross();
		pop += (*it)->getLocalPop();
	}
	cross -= getLocalCross();
	pop -= getLocalPop();

	Genotype child1(_p1->getGamete(), _p2->getGamete());
	gen++;
	cross++;
	pop += child1.computePop(nLoci, RM, gamma, _p1->getGenotype(), _p2->getGenotype());
	if (level == 0) 
		Node::printNode(out, child1, gen, pop, cross);

	if (_backcross == BackcrossP1)
	{
		Genotype child2(_p1->getGamete(), _gamete);
		gen++;
		cross++;
		pop += child2.computePop(nLoci, RM, gamma, _p1->getGenotype(), child1);

		if (level == 1)
			Node::printNode(out, child2, gen, pop, cross);

		if (_self)
		{
			Genotype child3(_gamete, _gamete);
			gen++;
			cross++;
			pop += child3.computePop(nLoci, RM, gamma, child2, child2);

			if (level == 2)
				Node::printNode(out, child3, gen, pop, cross);
		}
	}
	else if (_backcross == BackcrossP2)
	{
		Genotype child2(_p2->getGamete(), _gamete);
		gen++;
		cross++;
		pop += child2.computePop(nLoci, RM, gamma, _p2->getGenotype(), child1);

		if (level == 1)
			Node::printNode(out, child2, gen, pop, cross);

		if (_self)
		{
			Genotype child3(_gamete, _gamete);
			gen++;
			cross++;
			pop += child3.computePop(nLoci, RM, gamma, child2, child2);

			if (level == 2)
				Node::printNode(out, child3, gen, pop, cross);
		}
	}
	else if (_self)
	{
		Genotype child2(_gamete, _gamete);
		gen++;
		cross++;
		pop += child2.computePop(nLoci, RM, gamma, child1, child1);

		if (level == 1)
			Node::printNode(out, child2, gen, pop, cross);
	}
}
