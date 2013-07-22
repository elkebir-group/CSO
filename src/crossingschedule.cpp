/* 
 * crossingschedule.cpp
 *
 *  Created on: 28-apr-2011
 *      Author: M. El-Kebir
 */

#include "crossingschedule.h"
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include <algorithm>

CrossingSchedule::CrossingSchedule(const Data* pData)
  : _pData(pData)
  , _G()
  , _genotype(_G)
  , _ancestorSet(_G)
  , _pop(_G)
  , _gen(_G)
  , _cumPop(_G)
  , _cumCross(_G)
  , _upperChr(_G)
  , _targetNode(lemon::INVALID)
{
}

CrossingSchedule::CrossingSchedule(const CrossingSchedule& cs)
  : _pData(cs._pData)
  , _G()
  , _genotype(_G)
  , _ancestorSet(_G)
  , _pop(_G)
  , _gen(_G)
  , _cumPop(_G)
  , _cumCross(_G)
  , _upperChr(_G)
  , _targetNode(lemon::INVALID)
{
  DAG::NodeMap<Node> nr(cs._G);

  lemon::DigraphCopy<DAG, DAG> copier(cs._G, _G);
  copier.nodeMap(cs._genotype, _genotype);
  copier.nodeMap(cs._pop, _pop);
  copier.nodeMap(cs._gen, _gen);
  copier.nodeMap(cs._cumPop, _cumPop);
  copier.nodeMap(cs._cumCross, _cumCross);
  copier.arcMap(cs._upperChr, _upperChr);
  copier.node(cs._targetNode, _targetNode);
  copier.nodeRef(nr);
  copier.run();

  // build ancestor set
  for (NodeIt v(cs._G); v != lemon::INVALID; ++v)
  {
    const NodeSet& srcAnc = cs._ancestorSet[v];
    NodeSet& dstAnc = _ancestorSet[nr[v]];
    for (NodeSet::const_iterator it = srcAnc.begin(); it != srcAnc.end(); it++)
    {
      dstAnc.insert(nr[*it]);
    }
  }
}

void CrossingSchedule::updateAttributesDAG(Node node, BoolNodeMap& visited, bool updateGen)
{
  // we do a DFS here
  Node nodeP1 = lemon::INVALID;
  Node nodeP2 = lemon::INVALID;

  NodeSet& ancestors = _ancestorSet[node];
  ancestors.clear();

  for (InArcIt a(_G, node); a != lemon::INVALID; ++a)
  {
    if (_upperChr[a])
      nodeP1 = _G.source(a);
    else
      nodeP2 = _G.source(a);
  }

  if (nodeP1 != lemon::INVALID && !visited[nodeP1])
    updateAttributesDAG(nodeP1, visited, updateGen);
  if (nodeP2 != lemon::INVALID && !visited[nodeP2])
    updateAttributesDAG(nodeP2, visited, updateGen);

  if (nodeP1 == lemon::INVALID && nodeP2 == lemon::INVALID)
  {
    // node is a leaf
    _cumPop[node] = 0;
    _cumCross[node] = 0;
    if (updateGen)
      _gen[node] = 0;
  }
  else
  {
    NodeSet ancestorsP1 = nodeP1 != lemon::INVALID ? _ancestorSet[nodeP1] : NodeSet();
    NodeSet ancestorsP2 = nodeP2 != lemon::INVALID ? _ancestorSet[nodeP2] : NodeSet();
   
    // build ancestor set
    std::set_union(ancestorsP1.begin(), ancestorsP1.end(), 
      ancestorsP2.begin(), ancestorsP2.end(), 
      std::inserter(ancestors, ancestors.begin()));
    ancestors.insert(nodeP1);
    ancestors.insert(nodeP2);

    // update cumPop and cumCross
    _cumPop[node] = _pop[node];
    _cumCross[node] = 1;
    for (NodeSet::const_iterator it = ancestors.begin(); 
      it != ancestors.end(); it++)
    {
      _cumPop[node] += _pop[*it];
      if (!isLeaf(*it))
        _cumCross[node]++;
    }

    if (updateGen)
      _gen[node] = 1 + std::max(_gen[nodeP1], _gen[nodeP2]);
  }

  // flag as visited
  visited[node] = true;
}

void CrossingSchedule::recomputePop(Node node)
{
  if (isLeaf(node))
  {
    _pop[node] = 0;
  }
  else
  {
    Node nodeP1 = lemon::INVALID;
    Node nodeP2 = lemon::INVALID;

    for (InArcIt a(_G, node); a != lemon::INVALID; ++a)
    {
      if (_upperChr[a])
        nodeP1 = _G.source(a);
      else
        nodeP2 = _G.source(a);
    }

    _pop[node] = _genotype[node].computePop(_pData->getNumberOfLoci(),
        _pData->getRM(), _pData->getGamma(), _genotype[nodeP1], _genotype[nodeP2]);

    recomputePop(nodeP1);
    recomputePop(nodeP2);
  }
}

void CrossingSchedule::recomputePop()
{
  recomputePop(_targetNode);
  BoolNodeMap visited(_G, false);
  updateAttributesDAG(_targetNode, visited, false);
}

void CrossingSchedule::printNode(Node node, std::ostream& out) const
{
  static const Genotype& ideotype = _pData->getIdeotype();
  static const int nLoci = _pData->getNumberOfLoci();

  const unsigned long cumPop = _cumPop[node];
  const unsigned long gen = _gen[node];
  const unsigned long cumCross = _cumCross[node];
  const Genotype& genotype = _genotype[node];

  out << "\t\t";
  out << _G.id(node);
  //genotype.printGenotype(nLoci, false, out, "");
  out << " [fontcolor=white,style=rounded,shape=box,";
  out << "label=";

  out << "<<TABLE BORDER=\"0\" CELLSPACING=\"0\" CELLPADDING=\"4\" CELLBORDER=\"0\">" << std::endl;
  out << "\t\t\t<TR>" << std::endl;

  char buf[33];
  toBitstring(genotype.getC0(), nLoci, buf);
  for (int i = 0; i < nLoci; i++)
  {
    out << "\t\t\t\t<TD BGCOLOR=\"";

    if (ideotype(0, nLoci - i - 1) == genotype(0, nLoci - i - 1) ||
        ideotype(1, nLoci - i - 1) == genotype(0, nLoci - i - 1))
      out << "blue";
    else 
      out << "red";

    out << "\">" << buf[i] << "</TD>" << std::endl;
  }

  out << "\t\t\t\t<TD WIDTH=\"40\" ROWSPAN=\"2\"><FONT COLOR=\"black\">";
  out << cumPop << "<BR/>" << gen << "/" << cumCross << "<BR/>" 
    << _pData->getCostCrossover() * cumPop 
      + _pData->getCostGen() * gen + _pData->getCostNode() * cumCross;
  out << "</FONT></TD>" << std::endl;

  out << "\t\t\t</TR>" << std::endl;

  out << "\t\t\t<TR>" << std::endl;

  toBitstring(genotype.getC1(), nLoci, buf);
  for (int i = 0; i < nLoci; i++)
  {
    out << "\t\t\t\t<TD BGCOLOR=\"";

    if (ideotype(1, nLoci - i - 1) == genotype(1, nLoci - i - 1) ||
        ideotype(0, nLoci - i - 1) == genotype(1, nLoci - i - 1))
      out << "blue";
    else 
      out << "red";

    out << "\">" << buf[i] << "</TD>" << std::endl;
  }

  out << "\t\t\t</TR>" << std::endl;
  out << "\t\t</TABLE>>" << "];" << std::endl;
}

void CrossingSchedule::printNodes(std::ostream& out) const
{
  const unsigned long targetGen = _gen[_targetNode];
  std::vector<NodeList> nodesPerGen(targetGen + 1);

  const NodeSet& ancestors = _ancestorSet[_targetNode];
  for (NodeSet::const_iterator it = ancestors.begin(); 
    it != ancestors.end(); it++)
  {
    nodesPerGen[_gen[*it]].push_back(*it);
  }
  nodesPerGen[targetGen].push_back(_targetNode);

  for (unsigned long g = 0; g <= targetGen; g++)
  {
    const NodeList& nodes = nodesPerGen[g];

    if (nodes.size())
    {
      out << "\t{\n\t\trank = same;" << std::endl;

      for (NodeList::const_iterator it = nodes.begin(); it != nodes.end(); it++)
        printNode(*it, out);

      out << "\t}" << std::endl;
    }
  }
}

void CrossingSchedule::printEdges(Node node, 
                                  BoolNodeMap& visited, 
                                  std::ostream& out) const
{
  static const int nLoci = _pData->getNumberOfLoci();

  if (!visited[node])
  {
    for (InArcIt a(_G, node); a != lemon::INVALID; ++a)
    {
      Node parent = _G.source(a);
      out << '\t';
      out << _G.id(parent);
      //_genotype[parent].printGenotype(nLoci, false, out, "");
      out << " -> ";
      out << _G.id(node);
      //_genotype[node].printGenotype(nLoci, false, out, "");
      out << " [label=\"" << _pop[node] << "\"];" << std::endl;
      printEdges(parent, visited, out);
    }
    visited[node] = true;
  }
}

void CrossingSchedule::printDAG(std::ostream& out) const
{
  out << "digraph G {" << std::endl;
  printNodes(out);
  BoolNodeMap visited(_G, false);
  printEdges(_targetNode, visited, out);
  out << "}" << std::endl;
}

void CrossingSchedule::loadDAG(std::istream& in)
{
  lemon::digraphReader(_G, in)
    .node("target", _targetNode)
    .nodeMap("genotype", _genotype)
    .nodeMap("pop", _pop)
    .nodeMap("gen", _gen)
    .nodeMap("cumCross", _cumCross)
    .arcMap("upperChr", _upperChr)
    .run();

  BoolNodeMap visited(_G, false);
  updateAttributesDAG(_targetNode, visited, false);
}

void CrossingSchedule::saveDAG(std::ostream& out) const
{
  lemon::digraphWriter(_G, out)
    .node("target", _targetNode)
    .nodeMap("genotype", _genotype)
    .nodeMap("pop", _pop)
    .nodeMap("gen", _gen)
    .nodeMap("cumCross", _cumCross)
    .arcMap("upperChr", _upperChr)
    .run();
}
