/*
 * gaindividual.h
 *
 *  Created on: 23-mar-2009
 *      Author: M. El-Kebir
 */

#ifndef GAINDIVIDUAL_H_
#define GAINDIVIDUAL_H_

#include "csoga.h"
#include "node.h"

class GaIndividual
{
private:
	NodePtr _root;
	NodePtrSet _nodeSet;
	unsigned long _pop;
	unsigned long _cross;
	unsigned long _gen;
	int _nBadAlleles;

public:
	unsigned long getPop() const;
	unsigned long getCross() const;
	unsigned long getGen() const;
	const Genotype& getGenotype() const;
	int compare(const GaIndividual& individual) const;
	bool operator <(const GaIndividual& individual) const;
	bool operator <=(const GaIndividual& individual) const;
	bool operator >(const GaIndividual& individual) const;
	bool operator >=(const GaIndividual& individual) const;
	bool operator ==(const GaIndividual& individual) const;
	bool operator !=(const GaIndividual& individual) const;
	GaIndividual();
	GaIndividual(NodePtr root);
	GaIndividual(const GaIndividual& individual);
	~GaIndividual();
	GaIndividual& operator =(const GaIndividual& individual);
	static void crossover(GaIndividual& individual1, GaIndividual& individual2);
	void mutate();
	double getCost() const;
	void printDAG(std::ostream& out) const;
	double getFitness() const;
	NodePtr getRandomNode();
	void updateNodeSet();
	int getNrOfBadAlleles() const;
};

inline unsigned long GaIndividual::getPop() const
{
	return _pop;
}

inline unsigned long GaIndividual::getCross() const
{
	return _cross;
}

inline unsigned long GaIndividual::getGen() const
{
	return _gen;
}

inline int GaIndividual::compare(const GaIndividual& individual) const
{
	double cost = getCost();
	double individualCost = individual.getCost();

	if (cost < individualCost)
		return -1;
	else if (cost > individualCost)
		return 1;
	else
		return 0;
}

inline bool GaIndividual::operator <(const GaIndividual& individual) const
{
	return compare(individual) < 0;
}

inline bool GaIndividual::operator <=(const GaIndividual& individual) const
{
	return compare(individual) <= 0;
}

inline bool GaIndividual::operator >(const GaIndividual& individual) const
{
	return compare(individual) > 0;
}

inline bool GaIndividual::operator >=(const GaIndividual& individual) const
{
	return compare(individual) >= 0;
}

inline bool GaIndividual::operator ==(const GaIndividual& individual) const
{
	return compare(individual) == 0;
}

inline bool GaIndividual::operator !=(const GaIndividual& individual) const
{
	return compare(individual) != 0;
}

inline int GaIndividual::getNrOfBadAlleles() const
{
	return _nBadAlleles;
}

#endif
