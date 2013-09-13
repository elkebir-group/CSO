/*
 * leafnode.cpp
 *
 *  Created on: 20-mar-2009
 *      Author: M. El-Kebir
 */

#include "leafnode.h"
#include "genotypetable.h"
#include "data.h"

LeafNode::LeafNode(const Genotype& genotype, bool self, int gamete, double recc)
	: Node(gamete)
	, _parentalGenotype(genotype)
	, _offspringGenotype(genotype)
	, _self(self)
	, _localPop(1)
	, _recc(recc)
{
	update(false);
}

LeafNode::LeafNode(const LeafNode& leafNode)
	: Node(leafNode)
	, _parentalGenotype(leafNode._parentalGenotype)
	, _offspringGenotype(leafNode._offspringGenotype)
	, _self(leafNode._self)
	, _localPop(leafNode._localPop)
	, _recc(leafNode._recc)
{
}

LeafNode::LeafNode(const Genotype& genotype, double recc)
	: Node()
	, _parentalGenotype(genotype)
	, _offspringGenotype(genotype)
	, _self(randDouble() < 0.5)
	, _localPop(1)
	, _recc(recc)
{
	update(true);
}

LeafNode::~LeafNode()
{
}

void LeafNode::getAllNodes(NodePtrSet& nodeSet) const
{
	if (nodeSet.find(NodePtr(_this)) == nodeSet.end())
	{
		nodeSet.insert(NodePtr(_this));
	}
}

NodePtr LeafNode::swapNodes(const NodePtr oldNode, const NodePtr newNode)
{
	if (NodePtr(_this) == oldNode)
	{
		return newNode;	
	}
	else
	{
		return NodePtr(_this);
	}
}

void LeafNode::update(bool pickGamete)
{
	int nLoci = Data::getInstance()->getNumberOfLoci();
	double gamma = Data::getInstance()->getGamma();
	const DoubleMatrix& RM = Data::getInstance()->getRM();
	//int popMax = Data::getInstance()->getPopMax();

	_localPop = 1;
	if (pickGamete)
	{
		const GenotypeGamete& genotypeGamete = GenotypeTable::getInstance()->getGenotype(_parentalGenotype);
		_gamete = genotypeGamete.getGameteProportionate()._c;
	}
	
	if (_self)// && !_parentalGenotype.isHomozygous())
	{
		_offspringGenotype = Genotype(_gamete, _gamete);
		
		unsigned long pop = _offspringGenotype.computePop(nLoci, RM, gamma, _parentalGenotype, _parentalGenotype);
		//if (pop > popMax)
		//{
		//	_pOffspringGenotype = NULL;
		//	_self = false;
		//}
		//else
		{
			_localPop += pop;
		}
	}
	//else
	//{
	//	_self = false;
	//}
}

NodePtr LeafNode::create(const Genotype& genotype, bool self, int gamete, double recc)
{
	LeafNode* pLeafNode = new LeafNode(genotype, self, gamete, recc); 
	NodePtr res(pLeafNode);
	pLeafNode->_this = res;	
	return res;
}

NodePtr LeafNode::create(const Genotype& genotype, double recc)
{
	LeafNode* pLeafNode = new LeafNode(genotype, recc); 
	NodePtr res(pLeafNode);
	pLeafNode->_this = res;	
	return res;
}

NodePtr LeafNode::clone(const LeafNode* pLeaf)
{
	LeafNode* pClone = new LeafNode(*pLeaf);
	NodePtr res(pClone);
	pClone->_this = res;
	return res;
}
