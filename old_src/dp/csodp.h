/*
 * csodp.h
 *
 *  Created on: 23-mar-2009
 *      Author: M. El-Kebir
 */

#ifndef CSODP_H_
#define CSODP_H_

#include "cso.h"

#define CSO_DP_VERSION "23032009"

// forward declarations
class Solver;
class DpTable;
class DpItem;
class HeuristicCorrectLoci;
class HeuristicCorrectAlleles;
class HeuristicIndividualCost;
class HeuristicPureLinkage;
class HeuristicLargestSubGroupSize;
class HeuristicSubGroupSizes;

typedef std::set<DpItem> DpItemSet;
typedef std::list<DpItem> DpItemList;
typedef std::vector<DpItem> DpItemVector;
typedef std::list<const DpItem*> DpItemPointerList;
typedef std::stack<const DpItem*> DpItemPointerStack;

typedef std::tr1::unordered_map<int, std::tr1::unordered_map<int, DpItem> > DpMatrix;

typedef enum 
{
  HeuristicCorrectLociType = 1,
  HeuristicCorrectAllelesType = 2,
  //HeuristicIndividualCostType = 3,
  HeuristicLinkageType = 3,
  HeuristicLinkageLociType = 4,
  HeuristicLinkageAllelesType = 5,
  //HeuristicCorrectLociBonusPerLocusType = 6,
  //HeuristicCorrectAllelesBonusPerLocusType = 7,
  //HeuristicLinkageLociBonusType = 8,
  //HeuristicLinkageAllelesBonusType = 9,
  HeuristicPureLinkageType = 10,
  HeuristicLargestSubGroupSizeType = 11,
  HeuristicSubGroupSizesType = 12
} HeuristicType;

#endif
