/*
 * dptable.h
 *
 *  Created on: 09-mar-2009
 *      Author: s030858
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
	mutable int _useCount;

public:
	DpTable(int nLoci);
	virtual ~DpTable();
	virtual bool isPresent(const Genotype& genotype) const = 0;
	virtual DpItem& getItem(const Genotype& genotype) = 0;
	virtual const DpItem& getItem(const Genotype& genotype) const = 0;
	virtual DpItem& updateItem(const DpItem& item) = 0;
	int getUseCount() const;
};

inline int DpTable::getUseCount() const
{
	return _useCount;
}

#endif /* DPTABLE_H_ */