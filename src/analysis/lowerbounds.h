/*
 * lowerbounds.h
 *
 *  Created on: 3-may-2011
 *      Author: M. El-Kebir
 */

#ifndef LOWERBOUNDS_H_
#define LOWERBOUNDS_H_

#include "data.h"
#include "linkageanalysis.h"

class LowerBound
{
private:
  const Data* _pData;
  const bool _verbose;
  unsigned long _popLB;
  unsigned long _crossLB;
  unsigned long _genLB;
  unsigned long _popMaxLB;

  void compute();
  unsigned long solveSetCover();

public:
  LowerBound(const Data* pData, bool verbose);
  ~LowerBound() {}
  unsigned long getPopLB() const { return _popLB; }
  unsigned long getCrossLB() const { return _crossLB; }
  unsigned long getGenLB() const { return _genLB; }
  unsigned long getPopMaxLB() const { return _popMaxLB; }
};

#endif /* LOWERBOUNDS_H_ */
