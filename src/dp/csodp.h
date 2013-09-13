/*
 * csodp.h
 *
 *  Created on: 23-mar-2009
 *      Author: s030858
 */

#ifndef CSODP_H_
#define CSODP_H_

#include "cso.h"

#define CSO_DP_VERSION "23032009"

// forward declarations
class Solver;
class DpTable;
class DpItem;

typedef std::set<DpItem> DpItemSet;
typedef std::list<DpItem> DpItemList;
typedef std::vector<DpItem> DpItemVector;
typedef std::list<const DpItem*> DpItemPointerList;
typedef std::stack<const DpItem*> DpItemPointerStack;

typedef std::tr1::unordered_map<int, std::tr1::unordered_map<int, DpItem> > DpMatrix;

#endif
