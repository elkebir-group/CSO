/*
 * dptablearray.h
 *
 *  Created on: 09-mar-2009
 *      Author: M. El-Kebir
 */

#ifndef DPTABLEARRAY_H_
#define DPTABLEARRAY_H_

#include "dptable.h"

class DpTableArray : public DpTable
{
private:
  DpItem** _M;
  int getIndex(const Genotype& genotype) const;

public:
  DpTableArray(int crossOverCount = -1, bool restrictGametes = false);
  virtual ~DpTableArray();
  virtual bool isPresent(const Genotype& genotype) const;
  virtual const DpItem& getItem(const Genotype& genotype) const;
  virtual DpItem& getItem(const Genotype& genotype);
  virtual DpItem& updateItem(const DpItem& item);
};

inline int DpTableArray::getIndex(const Genotype& genotype) const
{
  int i = genotype.getC0();
  int j = genotype.getC1();

  return (1 << _nLoci) * i - (i + 1) * i / 2 + j;
}

#endif
