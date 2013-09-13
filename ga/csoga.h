/*
 * csoga.h
 *
 *  Created on: 23-mar-2009
 *      Author: M. El-Kebir
 */

#ifndef CSODP_H_
#define CSODP_H_

#include "cso.h"
//#include <boost/shared_ptr.hpp>
//#include <boost/weak_ptr.hpp>

#ifdef _MSC_VER
#include <memory>
#else
#include <tr1/memory>
#endif

#define CSO_GA_VERSION "24032009"

// forward declarations
class Node;
class LeafNode;
class InnerNode;
class GaIndividual;

//typedef boost::shared_ptr<Node> NodePtr;
//typedef boost::weak_ptr<Node> NodeWeakPtr;

typedef std::tr1::shared_ptr<Node> NodePtr;
typedef std::tr1::weak_ptr<Node> NodeWeakPtr;

typedef std::list<NodePtr> NodePtrList;
typedef std::set<NodePtr> NodePtrSet;

typedef std::set<NodeWeakPtr> NodeWeakPtrSet;

//typedef boost::shared_ptr<GaIndividual> GaIndividualPtr;

typedef std::tr1::shared_ptr<GaIndividual> GaIndividualPtr;

typedef std::vector<GaIndividual> GaIndividualVector;
typedef std::set<GaIndividual> GaIndividualSet;

#endif
