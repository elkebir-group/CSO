/*
 * data.h
 *
 *  Created on: 17-feb-2009
 *      Author: s030858
 */

#ifndef DATA_H_
#define DATA_H_

#include "cso.h"
#include <vector>
#include <limits>
#include "genotype.h"

class Data
{
private:
  Data(int nLoci,
       unsigned long popMax,
       int genMax,
       double gamma,
       double overallGamma,
       double costCrossover,
       double costGen,
       double costNode,
       bool allowSelfing,
       int c0,
       int c1);
  int _nLoci;
  unsigned long _popMax;
  int _genMax;
  DoubleMatrix _RM;
  DoubleMatrix _logRM;
  GenotypeSet _parents;
  Genotype _ideotype;
  double _gamma;
  double _overallGamma;
  double _costCrossover;
  double _costGen;
  double _costNode;
  double _probLowerBound;
  bool _allowSelfing;
  static Data* _pData;

  void computeLogRM();

public:
  virtual ~Data();
  static DoubleMatrix generateRM(const DoubleVector& cmVector);
  static Data* create(const char* fileName, bool readParents = true, bool allowPopMax = true);

  const DoubleMatrix& getRM() const;
  const DoubleMatrix& getLogRM() const;
  int getNumberOfLoci() const;
  const GenotypeSet& getParents() const;
  double getGamma() const;
  double getOverallGamma() const;
  unsigned long getPopMax() const;
  void printParents() const;
  void printRM() const;
  void printLogRM() const;
  double getCostCrossover() const;
  double getCostGen() const;
  double getCostNode() const;
  double getProbLowerBound() const;
  const Genotype& getIdeotype() const;
  double getCost(double crossover, unsigned long gen, unsigned long node) const;
  void resetPopMax();
  bool getAllowSelfing() const;
  static const Data* getInstance();
  void updateGamma(int nCrossings);
  int getGenMax() const;
  bool isIdeotypeHomozygous() const;
};

inline bool Data::isIdeotypeHomozygous() const
{
  return _pData->_ideotype.isHomozygous();
}

inline int Data::getGenMax() const
{
  return _genMax;
}

inline void Data::updateGamma(int nCrossings)
{
  if (0 < _overallGamma && _overallGamma < 1)
  {
    _gamma = pow(_overallGamma, 1. / nCrossings);
  }
}

inline void Data::resetPopMax()
{
  _popMax = std::numeric_limits<unsigned long>::max();
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

inline const DoubleMatrix& Data::getLogRM() const
{
  return _logRM;
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

inline double Data::getCostCrossover() const
{
  return _costCrossover;
}

inline double Data::getCostGen() const
{
  return _costGen;
}

inline double Data::getCostNode() const
{
  return _costNode;
}

inline double Data::getProbLowerBound() const
{
  return _probLowerBound;
}

inline const Genotype& Data::getIdeotype() const
{
  return _ideotype;
}

inline double Data::getCost(double crossover, unsigned long gen, unsigned long node) const
{
  return crossover * _costCrossover + gen * _costGen + node * _costNode;
}

inline bool Data::getAllowSelfing() const
{
  return _allowSelfing;
}

inline double Data::getOverallGamma() const
{
  return _overallGamma;
}

#endif /* DATA_H_ */
