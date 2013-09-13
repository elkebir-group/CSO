/*
 * gaindividual.h
 *
 *  Created on: 23-mar-2009
 *      Author: M. El-Kebir
 */

#include "gaindividual.h"
#include "innernode.h"
#include "leafnode.h"
#include "data.h"
#include <utility>

GaIndividual::GaIndividual()
	:  _root()
	, _nodeSet()
	, _pop(0)
	, _cross(0)
	, _gen(0)
	, _nBadAlleles(0)
{
	const Genotype& ideotype = Data::getInstance()->getIdeotype();
	const GenotypeSet& parents = Data::getInstance()->getParents();

	NodePtr p1 = Node::createSplitCross(parents, ideotype.getC0());
	if (p1->getGenotype() == ideotype)
	{
		_root = p1;
	}
	else
	{
		NodePtr p2 = Node::createSplitCross(parents, ideotype.getC1());
		if (p2->getGenotype() == ideotype)
		{
			_root = p2;
		}
		else
		{
			_root = InnerNode::create(p1, p2);
		}
	}

	updateNodeSet();
}

GaIndividual::GaIndividual(NodePtr root)
	: _root(root)
{
	updateNodeSet();
}

GaIndividual::GaIndividual(const GaIndividual& individual)
{
	// copy constructor
	_root = Node::clone(individual._root);
	updateNodeSet();
}

GaIndividual::~GaIndividual()
{
}

GaIndividual& GaIndividual::operator =(const GaIndividual& individual)
{
	if (this == &individual)
	{
		// self assignment
		return *this;
	}

	_root = Node::clone(individual._root);
	updateNodeSet();

	return *this;
}

void GaIndividual::updateNodeSet()
{
	_nodeSet.clear();
	_root->getAllNodes(_nodeSet);

	_root->updateGen();
	_gen = _root->getGen();
	_cross = _pop = 0;
	for (NodePtrSet::const_iterator it = _nodeSet.begin();
		it != _nodeSet.end(); it++)
	{
		_cross += (*it)->getLocalCross();
		_pop += (*it)->getLocalPop();
	}

	/* update _nBadAlleles */
	const Genotype& ideotype = Data::getInstance()->getIdeotype();
	const Genotype& rootGenotype = _root->getGenotype();
	int nLoci = Data::getInstance()->getNumberOfLoci();
	
	_nBadAlleles = 0;
	for (int i = 0; i < nLoci; i++)
	{
		if (rootGenotype(0, i) != ideotype(0, i))
			_nBadAlleles++;
		if (rootGenotype(1, i) != ideotype(1, i))
			_nBadAlleles++;
	}
}

NodePtr GaIndividual::getRandomNode()
{
	int size = (int) _nodeSet.size();
	assert(size);

	int index = randInt(size);
	NodePtrSet::iterator it = _nodeSet.begin();
	for (int i = 0; i < index; i++) it++;

	return *it;
}

void GaIndividual::crossover(GaIndividual& individual1, GaIndividual& individual2)
{
	NodePtr node1 = individual1.getRandomNode();
	NodePtr node2 = individual2.getRandomNode();

	NodePtr cpyNode1 = Node::clone(node1);
	NodePtr cpyNode2 = Node::clone(node2);

	individual1._root = individual1._root->swapNodes(node1, cpyNode2);
	individual2._root = individual2._root->swapNodes(node2, cpyNode1);

	individual1.updateNodeSet();
	individual2.updateNodeSet();
}

void GaIndividual::mutate()
{
	NodePtr randomNode = getRandomNode();

	InnerNode* pInnerNode = dynamic_cast<InnerNode*>(randomNode.get());
	if (pInnerNode && randDouble() < 0.8)
	{
		NodePtr parent1 = pInnerNode->getP1(), parent2 = pInnerNode->getP2();
		NodePtr newNode;

		if (pInnerNode->getGamete() == parent1->getGamete())
		{
			// remove bloat
			newNode = parent1;
		}
		else if (pInnerNode->getGamete() == parent2->getGamete())
		{
			// remove bloat
			newNode = parent2;
		}
		else if (randDouble() < 0.5)
		{
			newNode = InnerNode::create(parent1, parent2);
		}
		else if (dynamic_cast<InnerNode*>(parent1.get()) && randDouble() < 0.5)
		{
			InnerNode* pInnerNodeParent1 = dynamic_cast<InnerNode*>(parent1.get());
			if (randDouble() < 0.5)
			{
				// perform a backcross involving a random ancestor and parent2
				// subsequently cross the result with parent1
				NodePtrSet nodeSet;
				pInnerNodeParent1->getAllNodes(nodeSet);
				
				// pick a random ancestor
				int index = randInt((int) nodeSet.size());
				NodePtrSet::iterator it = nodeSet.begin();
				for (int i = 0; i < index; i++) it++;

				newNode = InnerNode::create(parent1, InnerNode::create(parent2, *it));
			}
			else
			{
				// cross grandparent2 with parent2
				// subsequently cross the result with grandparent1
				newNode = InnerNode::create(pInnerNodeParent1->getP1(), 
					InnerNode::create(pInnerNodeParent1->getP2(), parent2));
			}
		}
		else if (dynamic_cast<InnerNode*>(parent2.get()) && randDouble() < 0.5)
		{
			InnerNode* pInnerNodeParent2 = dynamic_cast<InnerNode*>(parent2.get());
			if (randDouble() < 0.5)
			{
				// perform a backcross involving a random ancestor and parent1
				// subsequently cross the result with parent2
				NodePtrSet nodeSet;
				pInnerNodeParent2->getAllNodes(nodeSet);
				
				// pick a random ancestor
				int index = randInt((int) nodeSet.size());
				NodePtrSet::iterator it = nodeSet.begin();
				for (int i = 0; i < index; i++) it++;

				newNode = InnerNode::create(parent2, InnerNode::create(parent1, *it));
			}
			else
			{
				// cross grandparent1 with parent1
				// subsequently cross the result with grandparent2
				newNode = InnerNode::create(pInnerNodeParent2->getP2(), 
					InnerNode::create(pInnerNodeParent2->getP1(), parent1));
			}
		}
		else
		{
			newNode = InnerNode::create(parent1, parent2);
		}
		
		_root = _root->swapNodes(randomNode, newNode);
		updateNodeSet();
	}
}

const Genotype& GaIndividual::getGenotype() const
{
	return _root->getGenotype();
}


double GaIndividual::getCost() const
{
	return Data::getInstance()->getCost(_pop, _gen, _cross);
}

double GaIndividual::getFitness() const
{
	return getCost() + _nBadAlleles * 500000;
}

void GaIndividual::printDAG(std::ostream& out) const
{
	out << "digraph G {" << std::endl;

	std::vector<std::list<std::pair<NodePtr, int> > > nodesPerGeneration(_gen + _root->getLocalGen());

	for (NodePtrSet::const_iterator it = _nodeSet.begin();
		it != _nodeSet.end(); it++)
	{
		const NodePtr& node = *it;
		int localGen = (int) node->getLocalGen();
		for (int i = 0; i < localGen; i++)
		{
			nodesPerGeneration[node->getGen() + i].push_back(std::make_pair(node, i));
		}
	}

	for (unsigned long i = 0; i < _gen + _root->getLocalGen(); i++)
	{
		const std::list<std::pair<NodePtr, int> >& list = nodesPerGeneration[i];

		if (list.size() != 0)
		{
			out << "\t{\n\t\trank = same;\n";
		}

		for (std::list<std::pair<NodePtr, int> >::const_iterator it = list.begin(); it != list.end(); it++)
		{
			it->first->printNode(out, it->second);
		}

		if (list.size() != 0)
		{
			out << "\t}\n";
		}
	}

	for (NodePtrSet::const_iterator it = _nodeSet.begin();
		it != _nodeSet.end(); it++)
	{
		(*it)->printEdges(out);
	}

	out << "}" << std::endl;
}
