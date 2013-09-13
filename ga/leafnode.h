/*
 * leafnode.h
 *
 *  Created on: 20-mar-2009
 *      Author: M. El-Kebir
 */

#ifndef LEAFNODE_H_
#define LEAFNODE_H_

#include "node.h"

/* TODO
 * uitbreiding: 
 *   - nieuwe member _cross
 *   - nieuwe member _gen
 */

class LeafNode : public Node
{
private:
	const Genotype _parentalGenotype;
	Genotype _offspringGenotype;
	bool _self;
	unsigned long _localPop;
	double _recc;

	void update(bool pickGamete);
	LeafNode(const Genotype& genotype, bool self, int gamete, double recc);
	LeafNode(const LeafNode& leafNode);
	LeafNode(const Genotype& genotype, double recc);

public:
	static NodePtr create(const Genotype& genotype, bool self, int gamete, double recc);
	static NodePtr create(const Genotype& genotype, double recc);
	static NodePtr clone(const LeafNode* pLeaf);
	~LeafNode();
	void updateGen();
	unsigned long getGen() const;
	unsigned long getLocalGen() const;
	unsigned long getLocalCross() const;
	unsigned long getLocalPop() const;
	double getRecc() const;
	const Genotype& getGenotype() const;
	void getAllNodes(NodePtrSet& nodeSet) const;
	NodePtr swapNodes(const NodePtr oldNode, const NodePtr newNode);
	void printEdges(std::ostream& out) const;
	void printNode(std::ostream& out, int level) const;
};

inline const Genotype& LeafNode::getGenotype() const
{
	if (_self)
	{
		return _offspringGenotype;
	}
	else
	{
		return _parentalGenotype;
	}
}

inline unsigned long LeafNode::getLocalGen() const
{
	return 1 + (_self ? 1 : 0);
}

inline unsigned long LeafNode::getLocalCross() const
{
	return (_self ? 1 : 0);
}

inline unsigned long LeafNode::getLocalPop() const
{
	return _localPop;
}

inline unsigned long LeafNode::getGen() const
{
	return (_self ? 1 : 0);
}

inline double LeafNode::getRecc() const
{
	return _recc;
}

inline void LeafNode::updateGen()
{
}

inline void LeafNode::printEdges(std::ostream& out) const
{
	if (_self)
	{
		printEdge(out, _parentalGenotype, _parentalGenotype, _offspringGenotype);
	}
}

inline void LeafNode::printNode(std::ostream& out, int level) const
{
	assert(0 <= level && level < (int) getLocalGen());

	if (level == 0)
	{
		Node::printNode(out, _parentalGenotype, 0, 1, 0);
	}
	else if (level == 1)
	{
		Node::printNode(out, _offspringGenotype, 1, _localPop, 1);
	}
}

#endif
