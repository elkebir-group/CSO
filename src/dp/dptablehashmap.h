/*
 * dptablehashmap.h
 *
 *  Created on: 09-mar-2009
 *      Author: s030858
 */

#ifndef DPTABLEHASHMAP_H_
#define DPTABLEHASHMAP_H_

#include "dptable.h"

class DpTableHashMap : public DpTable
{
private:
	mutable DpMatrix _M;

public:
	DpTableHashMap(int nLoci);
	virtual ~DpTableHashMap();
	virtual bool isPresent(const Genotype& genotype) const;
	virtual const DpItem& getItem(const Genotype& genotype) const;
	virtual DpItem& getItem(const Genotype& genotype);
	virtual DpItem& updateItem(const DpItem& item);
};

#endif