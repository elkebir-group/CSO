/* 
 * crossingschedule.h
 *
 *  Created on: 28-apr-2011
 *      Author: M. El-Kebir
 */

#ifndef CROSSINGSCHEDULE_H_
#define CROSSINGSCHEDULE_H_

#include <lemon/smart_graph.h>
#include <set>
#include <vector>
#include <list>
#include <ostream>
#include <assert.h>
#include "common/data.h"
#include "common/genotype.h"

class CrossingSchedule
{
public:
  /// Crossing schedule type
  typedef lemon::SmartDigraph DAG;

protected:
  DIGRAPH_TYPEDEFS(DAG);
  typedef std::set<Node> NodeSet;
  typedef std::list<Node> NodeList;
  typedef DAG::ArcMap<double> ProbMap;
  typedef DAG::ArcMap<bool> UpperChrMap;
  typedef DAG::NodeMap<NodeSet> AncestorMap;

  typedef DAG::NodeMap<Genotype> GenotypeMap;
  typedef DAG::NodeMap<double> PopMap;
  typedef DAG::NodeMap<unsigned long> GenMap;
  typedef DAG::NodeMap<double> CumPopMap;
  typedef DAG::NodeMap<unsigned long> CumCrossMap;
  typedef DAG::NodeMap<int> IdxMap;

  const Data* _pData;
  DAG _G;
  GenotypeMap _genotype;
  AncestorMap _ancestorSet;
  PopMap _pop;
  PopMap _prob;
  GenMap _gen;
  CumPopMap _cumPop;
  CumCrossMap _cumCross;
  UpperChrMap _upperChr;
  Node _targetNode;

  void updateAttributesDAG(Node node, BoolNodeMap& visited, bool updateGen);
  void recomputePop(Node node);
  bool isLeaf(Node node) const;
  void printNode(Node node, std::ostream& out) const;
  void printNodes(std::ostream& out) const;
  virtual void printEdges(Node node, BoolNodeMap& visited, std::ostream& out, bool prob) const;

public:
  CrossingSchedule(const Data* pData);
  CrossingSchedule(const CrossingSchedule& cs);
  virtual ~CrossingSchedule() {}
  void printDAG(std::ostream& out, bool prob = false) const;
  bool isFeasible() const;
  double getPop() const;
  unsigned long getGen() const;
  unsigned long getCross() const;
  double getCost() const;
  void loadDAG(std::istream& in);
  void saveDAG(std::ostream& out) const;
  void recomputePop();
};

inline bool CrossingSchedule::isFeasible() const
{
  return _targetNode != lemon::INVALID;
}

inline double CrossingSchedule::getPop() const
{
  return isFeasible() ? _cumPop[_targetNode] : -1;
}

inline unsigned long CrossingSchedule::getCross() const
{
  return isFeasible() ? _cumCross[_targetNode] : -1;
}

inline unsigned long CrossingSchedule::getGen() const
{
  return isFeasible() ? _gen[_targetNode] : -1;
}

inline double CrossingSchedule::getCost() const
{
  return _pData->getCostCrossover() * getPop()
    + _pData->getCostGen() *  getGen()
    + _pData->getCostNode() * getCross();
}

inline bool CrossingSchedule::isLeaf(Node node) const
{
  return InArcIt(_G, node) == lemon::INVALID;
}

#endif /* CROSSINGSCHEDULE_H_ */
