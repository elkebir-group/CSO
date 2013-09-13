/*
 * node.h
 *
 *  Created on: 20-mar-2009
 *      Author: M. El-Kebir
 */

#ifndef NODE_H_
#define NODE_H_

#include "cso.h"
#include "csoga.h"
#include "genotype.h"
#include "gaindividual.h"

class Node
{
protected:
	int _gamete;
	/* nodig voor swapNodes */
	NodeWeakPtr _this;

protected:
	Node();
	Node(int gamete);
	virtual ~Node();

public:
	static NodePtr clone(const Node* pNode);
	static NodePtr clone(const NodePtr pNode);
	static NodePtr createSplitCross(const GenotypeSet& parents, int targetGamete);
	virtual void updateGen() = 0;
	virtual unsigned long getGen() const = 0;
	virtual unsigned long getLocalGen() const = 0;
	virtual unsigned long getLocalCross() const = 0;
	virtual unsigned long getLocalPop() const = 0;
	virtual double getRecc() const = 0;
	int getGamete() const;
	virtual const Genotype& getGenotype() const = 0;
	virtual void getAllNodes(NodePtrSet& nodeSet) const = 0;
	virtual NodePtr swapNodes(const NodePtr oldNode, const NodePtr newNode) = 0;
	virtual void printEdges(std::ostream& out) const = 0;
	virtual void printNode(std::ostream& out, int level) const = 0;
	static void printNode(std::ostream& out, const Genotype& genotype, unsigned long gen, unsigned long pop, unsigned long cross);
	static void printEdge(std::ostream& out, const Genotype& parent1, const Genotype& parent2, const Genotype& child);
};

inline Node::Node()
	: _gamete(INVALID_GAMETE)
	, _this()
{
}

inline Node::Node(int gamete)
	: _gamete(gamete)
	, _this()
{
}

inline Node::~Node()
{
}

inline int Node::getGamete() const
{
	return _gamete;
}

#endif
