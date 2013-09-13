/*
 * dptable.h
 *
 *  Created on: 09-mar-2009
 *      Author: M. El-Kebir
 */

#ifndef DPTABLE_H_
#define DPTABLE_H_

#include "csodp.h"
#include "dpitem.h"
#include "genotype.h"

class DpTable
{
protected:
	const int _nLoci;
	const int _nGenotypes;
	const DoubleMatrix _RM;
	const double _probLowerBound;
	const Genotype _ideotype;
	const int _crossOverCount;
	const bool _restrictGametes;
	mutable int _useCount;

public:
	DpTable(int crossOverCount, bool restrictGametes);
	virtual ~DpTable();
	virtual bool isPresent(const Genotype& genotype) const = 0;
	virtual DpItem& getItem(const Genotype& genotype) = 0;
	virtual const DpItem& getItem(const Genotype& genotype) const = 0;
	virtual DpItem& updateItem(const DpItem& item) = 0;
	int getUseCount() const;
	int getCrossOverCount() const;
	bool getRestrictGametes() const;
};

inline int DpTable::getUseCount() const
{
	return _useCount;
}

inline int DpTable::getCrossOverCount() const
{
	return _crossOverCount;
}

inline bool DpTable::getRestrictGametes() const
{
	return _restrictGametes;
}

#endif /* DPTABLE_H_ */