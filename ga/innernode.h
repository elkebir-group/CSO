/*
 * innernode.h
 *
 *  Created on: 20-mar-2009
 *      Author: M. El-Kebir
 */

#ifndef INNERNODE_H_
#define INNERNODE_H_

#include "node.h"

typedef enum
{
	NoBackcross = 0,
	BackcrossP1 = 1,
	BackcrossP2 = 2
} BackcrossType;

class InnerNode : public Node
{
private:
	Genotype _genotype;
	NodePtr _p1;
	NodePtr _p2;
	bool _self;
	BackcrossType _backcross;
	unsigned long _localPop;
	unsigned long _gen;

	void update(bool pickGamete);
	InnerNode(const InnerNode& innerNode);
	InnerNode(NodePtr p1, NodePtr p2);
	InnerNode(NodePtr p1, NodePtr p2, bool self, BackcrossType backcross, int gamete);

public:
	static NodePtr create(NodePtr p1, NodePtr p2);
	static NodePtr create(NodePtr p1, NodePtr p2, bool self, BackcrossType backcross, int gamete);
	static NodePtr clone(const InnerNode* pInnerNode);
	~InnerNode();
	void updateGen();
	unsigned long getGen() const;
	unsigned long getLocalGen() const;
	unsigned long getLocalCross() const;
	unsigned long getLocalPop() const;
	double getRecc() const;
	const Genotype& getGenotype() const;
	InnerNode& operator =(const InnerNode& innerNode);
	void getAllNodes(NodePtrSet& nodeSet) const;
	NodePtr swapNodes(const NodePtr oldNode, const NodePtr newNode);
	void printEdges(std::ostream& out) const;
	void printNode(std::ostream& out, int level) const;
	NodePtr getP1();
	NodePtr getP2();
};

inline NodePtr InnerNode::getP1()
{
	return _p1;
}

inline NodePtr InnerNode::getP2()
{
	return _p2;
}

inline const Genotype& InnerNode::getGenotype() const
{
	return _genotype;
}

inline unsigned long InnerNode::getLocalCross() const
{
	return 1 + (_self ? 1 : 0) + (_backcross != NoBackcross ? 1 : 0);
}

inline unsigned long InnerNode::getLocalGen() const
{
	return 1 + (_self ? 1 : 0) + (_backcross != NoBackcross ? 1 : 0);
}

inline unsigned long InnerNode::getLocalPop() const
{
	return _localPop;
}

inline unsigned long InnerNode::getGen() const
{
	return _gen;
}

#endif
