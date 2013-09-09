/*
 * data.h
 *
 *  Created on: 17-feb-2009
 *      Author: M. El-Kebir
 */

#ifndef DATA_H_
#define DATA_H_

#include "cso.h"
#include <vector>
#include "genotype.h"

class Data
{
private:
	Data(int nLoci, unsigned long popMax, double gamma, double costPop, double costGen, double costCross, bool allowSelfing, int c0, int c1);
	int _nLoci;
	unsigned long _popMax;
	DoubleMatrix _RM;
	GenotypeSet _parents;
	Genotype _ideotype;
	double _gamma;
	double _costPop;
	double _costGen;
	double _costCross;
	double _probLowerBound;
	bool _allowSelfing;
	static Data* _pData;

public:
	virtual ~Data();
	static DoubleMatrix generateRM(const DoubleVector& cmVector);
	static Data* create(const char* fileName, bool readParents = true, bool allowPopMax = true);
	//static Data* create(const int nLoci, const GenotypeSet& parents, const Genotype& ideotype, 
	//	double gamma, double costPop, double costGen, double costCross, unsigned long popMax, bool allowSelfing);

	const DoubleMatrix& getRM() const;
	int getNumberOfLoci() const;
	const GenotypeSet& getParents() const;
	double getGamma() const;
	unsigned long getPopMax() const;
	void printParents() const;
	void printRM() const;
	double getCostPop() const;
	double getCostGen() const;
	double getCostCross() const;
	double getProbLowerBound() const;
	const Genotype& getIdeotype() const;
	double getCost(unsigned long pop, unsigned long gen, unsigned long cross) const;
	void resetPopMax();
	bool getAllowSelfing() const;
	static const Data* getInstance();
};

inline void Data::resetPopMax()
{
	_popMax = ULONG_MAX;
	_probLowerBound = 0;
}

inline const Data* Data::getInstance()
{
	assert(_pData);
	return _pData;
}

inline const DoubleMatrix& Data::getRM() const
{
	return _RM;
}

inline int Data::getNumberOfLoci() const
{
	return _nLoci;
}

inline const GenotypeSet& Data::getParents() const
{
	return _parents;
}

inline double Data::getGamma() const
{
	return _gamma;
}

inline unsigned long Data::getPopMax() const
{
	return _popMax;
}

inline double Data::getCostPop() const
{
	return _costPop;
}

inline double Data::getCostGen() const
{
	return _costGen;
}

inline double Data::getCostCross() const
{
	return _costCross;
}

inline double Data::getProbLowerBound() const
{
	return _probLowerBound;
}

inline const Genotype& Data::getIdeotype() const
{
	return _ideotype;
}

inline double Data::getCost(unsigned long pop, unsigned long gen, unsigned long cross) const
{
	return pop * _costPop + gen * _costGen + cross * _costCross;
}

inline bool Data::getAllowSelfing() const
{
	return _allowSelfing;
}

#endif /* DATA_H_ */
