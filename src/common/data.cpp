/*
 * data.cpp
 *
 *  Created on: 17-feb-2009
 *      Author: s030858
 */

#include "data.h"
#include "genotypegamete.h"
#include "tinyxml/tinyxml.h"
#include <iostream>
#include <string>
#include <sstream>
#include <limits.h>
#include <limits>
#include <math.h>

Data* Data::_pData = NULL;

Data::Data(int nLoci,
           unsigned long popMax,
           int genMax,
           double gamma,
           double overallGamma,
           double costCrossover,
           double costGen,
           double costNode,
           bool allowSelfing,
           int c0,
           int c1)
  : _nLoci(nLoci)
  , _popMax(popMax)
  , _genMax(genMax)
  , _RM(nLoci)
  , _logRM(nLoci)
  , _ideotype(c0, c1)
  , _gamma(gamma)
  , _overallGamma(overallGamma)
  , _costCrossover(costCrossover)
  , _costGen(costGen)
  , _costNode(costNode)
  , _probLowerBound(popMax == std::numeric_limits<unsigned long>::max() ? 0 : 1.0 - pow((1.0 - _gamma), 1.0 / _popMax))
  , _allowSelfing(allowSelfing)
{
  for (int i = 0; i < nLoci; i++)
  {
    _RM[i] = std::vector<double>(nLoci);
    _logRM[i] = std::vector<double>(nLoci);
  }
}

Data::~Data()
{
  _pData = NULL;
}

DoubleMatrix Data::generateRM(const DoubleVector& cmVector)
{
  const int n = (int) cmVector.size();
  
  DoubleMatrix res(n);
  for (int i = 0; i < n; i++)
    res[i] = DoubleVector(n);
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (i < j)
      {
        res[i][j] = 0.5 * (1 - exp(-2.0 * (cmVector[j] - cmVector[i]) / 100.0));
      }
      else if (i > j)
      {
        res[i][j] = 0.5 * (1 - exp(-2.0 * (cmVector[i] - cmVector[j]) / 100.0));
      }
      else
      {
        res[i][j] = 0;
      }
    }
  }

  return res;
}

void Data::computeLogRM()
{
  for (int i = 0; i < _nLoci; i++)
  {
    for (int j = 0; j < _nLoci; j++)
    {
      _logRM[i][j] = log(_RM[i][j]);
    }
  }
}

Data* Data::create(const char* fileName, bool readParents, bool allowPopMax)
{
  TiXmlDocument xmlDoc(fileName);
  if (!xmlDoc.LoadFile())
  {
    std::cerr << xmlDoc.ErrorDesc() << std::endl;
    return NULL;
  }

  TiXmlElement* pRootElement = xmlDoc.RootElement();
  if (!pRootElement || strcmp(pRootElement->Value(), "CSO"))
  {
    std::cerr << "Expected root element with name 'CSO'" << std::endl;
    return NULL;
  }

  /* Get constants */
  double costCrossover = 1, costGen = 1, costNode = 1;
  pRootElement->Attribute("costCrossover", &costCrossover);
  pRootElement->Attribute("costGen", &costGen);
  pRootElement->Attribute("costNode", &costNode);

  /* Get number of loci */
  int nLoci;
  if (!pRootElement->Attribute("nLoci", &nLoci) || nLoci < 0)
  {
    std::cerr << "'/CSO/@nLoci' should be a positive integer" << std::endl;
    return NULL;
  }

  /* Get gamma or overall gamma */
  double gamma = -1, overallGamma = -1;
  if (!(pRootElement->Attribute("gamma", &gamma) && (0 < gamma && gamma < 1))
      && !(pRootElement->Attribute("overallGamma", &overallGamma) && (0 < overallGamma && overallGamma < 1)))
  {
    std::cerr << "Either '/CSO/@overallGamma' or '/CSO/@gamma' need to be specified\n"
              << "and they need to be in the range (0,1)" << std::endl;
    return NULL;
  }

  /* Get genMax */
  int genMax = std::numeric_limits<int>::max();
  if (pRootElement->Attribute("genMax", &genMax) && genMax < 0)
  {
    std::cerr << "'/CSO/@genMax' should be a positive integer" << std::endl;
    return NULL;
  }

  /* Get selfing */
  int selfing = 1;
  if (pRootElement->Attribute("selfing", &selfing) && selfing != 0 && selfing != 1)
  {
    std::cerr << "'/CSO/@selfing' should be 0 or 1" << std::endl;
    return NULL;
  }

  unsigned long popMax = std::numeric_limits<unsigned long>::max();
  if (allowPopMax)
  {
    int tempPopMax = std::numeric_limits<int>::max();
    pRootElement->Attribute("popMax", &tempPopMax);
    if (tempPopMax < 0)
    {
      std::cerr << "'/CSO/@popMax' should be in the range [0,\\infty)" << std::endl;
      return NULL;
    }
    else if (tempPopMax != std::numeric_limits<int>::max())
    {
      popMax = tempPopMax;
    }
  }

  /* Get ideotype */
  int c0, c1;
  TiXmlElement* pIdeotypeElement = pRootElement->FirstChildElement("Ideotype");
  if (!pIdeotypeElement)
  {
    std::cerr << "Missing mandatory element '/CSO/Ideotype'" << std::endl;
    return NULL;
  }

  const char* buf;
  if (!(buf = pIdeotypeElement->Attribute("c0")))
  {
    std::cerr << "Missing mandatory attribute '/CSO/Ideotype/@c0'" << std::endl;
    return NULL;
  }
  c0 = fromBitstring(nLoci, buf);

  if (!(buf = pIdeotypeElement->Attribute("c1")))
  {
    std::cerr << "Missing mandatory attribute '/CSO/Ideotype/@c1'" << std::endl;
    return NULL;
  }
  c1 = fromBitstring(nLoci, buf);

  Data::_pData = new Data(nLoci,
                          popMax,
                          genMax,
                          gamma,
                          overallGamma,
                          costCrossover,
                          costGen,
                          costNode,
                          selfing != 0,
                          c0,
                          c1);

  /* Get RM */
  if (pRootElement->FirstChildElement("cM"))
  {
    TiXmlElement* pCMElement = pRootElement->FirstChildElement("cM");

    std::string string(pCMElement->GetText());
    std::istringstream stream(string);

    std::vector<double> cmVector(nLoci);
    for (int i = 0; i < nLoci; i++)
    {
      stream >> cmVector[i];
    }

    Data::_pData->_RM = generateRM(cmVector);
  }
  else
  {
    int row = 0;
    for (TiXmlElement* pRMElement = pRootElement->FirstChildElement("RM");
      pRMElement;
      pRMElement = pRMElement->NextSiblingElement())
    {
      std::string string(pRMElement->GetText());
      std::istringstream stream(string);

      for (int i = 0; i < nLoci; i++)
      {
        stream >> Data::_pData->_RM[row][i];
      }

      row++;
    }
  }

  Data::_pData->computeLogRM();

  /* Get parents */
  if (readParents)
  {
    TiXmlElement* pParentsElement = pRootElement->FirstChildElement("Parents");
    if (!pParentsElement)
    {
      std::cerr << "Missing mandatory element '/CSO/Parents'" << std::endl;
      delete Data::_pData;
      Data::_pData = NULL;
      return NULL;
    }

    for (TiXmlElement* pParentElement = pParentsElement->FirstChildElement("Parent");
      pParentElement;
      pParentElement = pParentElement->NextSiblingElement())
    {
      int c0, c1;

      if (!(buf = pParentElement->Attribute("c0")))
      {
        std::cerr << "Missing mandatory attribute '/CSO/Parents/Parent/@c0'" << std::endl;
        delete Data::_pData;
        Data::_pData = NULL;
        return NULL;
      }
      c0 = fromBitstring(nLoci, buf);

      if (!(buf = pParentElement->Attribute("c1")))
      {
        std::cerr << "Missing mandatory attribute '/CSO/Parents/Parent/@c1'" << std::endl;
        delete Data::_pData;
        Data::_pData = NULL;
        return NULL;
      }
      c1 = fromBitstring(nLoci, buf);

      GenotypeGamete parent(c0, c1);
      //parent.computeGametes(nLoci, Data::_pData->_RM, Data::_pData->_probLowerBound);
      //parent.computeGametesCumulative();
      //if (parent.getGametes().size())
      //{
      Data::_pData->_parents.insert(parent);
      //}
      //else
      //{
      //  std::cerr << "Skipped parent ";
      //  parent.printGenotype(nLoci, false, std::cerr);
      //  std::cerr << ", as it produces no gametes within N_{max}." << std::endl;
      //}
    }
  }

  return Data::_pData;
}

void Data::printParents() const
{
  int n = (int) _parents.size();

  std::cout << "Number of parents: " << n << std::endl;

  for (GenotypeSet::const_iterator it = _parents.begin(); it != _parents.end(); it++)
  {
    it->printGenotype(_nLoci);
  }
}

void Data::printRM() const
{
  std::cout << "// Number of loci: " << _nLoci << std::endl;

  for (int i = 0; i < _nLoci; i++)
  {
    std::cout << "// ";
    for (int j = 0; j < _nLoci; j++)
    {
      std::cout << _RM[i][j] << "\t";
    }
    std::cout << std::endl;
  }
}

void Data::printLogRM() const
{
  std::cout << "// Number of loci: " << _nLoci << std::endl;

  for (int i = 0; i < _nLoci; i++)
  {
    std::cout << "// ";
    for (int j = 0; j < _nLoci; j++)
    {
      std::cout << _logRM[i][j] << "\t";
    }
    std::cout << std::endl;
  }
}
