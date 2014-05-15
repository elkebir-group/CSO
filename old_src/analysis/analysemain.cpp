/*
 * analysemain.cpp
 *
 *  Created on: 14-apr-2009
 *      Author: M. El-Kebir
 */

#define CSO_ANALYSIS_VER "14042009"

#include <stdio.h>
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <exception>
#include <time.h>
#include "genotype.h"
#include "data.h"
#include "cso.h"
#include "linkageanalysis.h"
#include "coveranalysis.h"

namespace po = boost::program_options;

int analyse(const char* pInputFileName)
{
  Data* pData = Data::create(pInputFileName);
  if (!pData) return 1;

  int ideotype0 = pData->getIdeotype().getC0();
  int ideotype1 = pData->getIdeotype().getC1();

  LinkageAnalysis linkageAnalysis(pData);
  CoverAnalysis coverAnalysis(pData);

  const GenotypePopVector& popVectorIdeotype0 = coverAnalysis.getGenotypePopVector(ideotype0);
  for (int i = 0; i < pData->getNumberOfLoci(); i++)
  {
    std::cout << "Locus " << i << ": ";
    popVectorIdeotype0[i]._genotype.printGenotype(pData->getNumberOfLoci(), false);
    std::cout << " - " << popVectorIdeotype0[i]._pop << " / " << popVectorIdeotype0[i]._count << std::endl;
  }

  std::cout << "Chromosome 0: unlinked alleles" << std::endl;
  const LociPairList& lociPairDifferentC0 = linkageAnalysis.getUnlinkedLoci(ideotype0);
  for (LociPairList::const_iterator it = lociPairDifferentC0.begin(); it != lociPairDifferentC0.end(); it++)
  {
    LinkagePopPair linkagePopPair = linkageAnalysis.getLinkagePopPair(ideotype0, *it);
    printf("(%d,%d)\t%lu\n", it->first, it->second, linkagePopPair.second);
  }

  std::cout << "Chromosome 0: weakly linked alleles" << std::endl;
  const LociPairList& lociPairSameC0 = linkageAnalysis.getWeaklyLinkedLoci(ideotype0);
  for (LociPairList::const_iterator it = lociPairSameC0.begin(); it != lociPairSameC0.end(); it++)
  {
    LinkagePopPair linkagePopPair = linkageAnalysis.getLinkagePopPair(ideotype0, *it);
    printf("(%d,%d)\t%lu\n", it->first, it->second, linkagePopPair.second);
  }

  std::cout << "Chromosome 0: linked alleles" << std::endl;
  const LociPairList& lociPairLinkedC0 = linkageAnalysis.getLinkedLoci(ideotype0);
  for (LociPairList::const_iterator it = lociPairLinkedC0.begin(); it != lociPairLinkedC0.end(); it++)
  {
    LinkagePopPair linkagePopPair = linkageAnalysis.getLinkagePopPair(ideotype0, *it);
    printf("(%d,%d)\t%lu\n", it->first, it->second, linkagePopPair.second);
  }

  if (pData->getIdeotype().isHomozygous())
    return 0;

  const GenotypePopVector& popVectorIdeotype1 = coverAnalysis.getGenotypePopVector(ideotype1);
  for (int i = 0; i < pData->getNumberOfLoci(); i++)
  {
    std::cout << "Locus " << i << ": ";
    popVectorIdeotype1[i]._genotype.printGenotype(pData->getNumberOfLoci(), false);
    std::cout << " - " << popVectorIdeotype1[i]._pop << " / " << popVectorIdeotype0[i]._count << std::endl;
  }

  std::cout << "Chromosome 1: unlinked loci in different parents" << std::endl;
  const LociPairList& lociPairDifferentC1 = linkageAnalysis.getUnlinkedLoci(ideotype1);
  for (LociPairList::const_iterator it = lociPairDifferentC1.begin(); it != lociPairDifferentC1.end(); it++)
  {
    LinkagePopPair linkagePopPair = linkageAnalysis.getLinkagePopPair(ideotype1, *it);
    printf("(%d,%d)\t%lu\n", it->first, it->second, linkagePopPair.second);
  }

  std::cout << "Chromosome 1: unlinked loci in same parent" << std::endl;
  const LociPairList& lociPairSameC1 = linkageAnalysis.getWeaklyLinkedLoci(ideotype1);
  for (LociPairList::const_iterator it = lociPairSameC1.begin(); it != lociPairSameC1.end(); it++)
  {
    LinkagePopPair linkagePopPair = linkageAnalysis.getLinkagePopPair(ideotype1, *it);
    printf("(%d,%d)\t%lu\n", it->first, it->second, linkagePopPair.second);
  }

  std::cout << "Chromosome 1: linked alleles" << std::endl;
  const LociPairList& lociPairLinkedC1 = linkageAnalysis.getLinkedLoci(ideotype1);
  for (LociPairList::const_iterator it = lociPairLinkedC1.begin(); it != lociPairLinkedC1.end(); it++)
  {
    LinkagePopPair linkagePopPair = linkageAnalysis.getLinkagePopPair(ideotype1, *it);
    printf("(%d,%d)\t%lu\n", it->first, it->second, linkagePopPair.second);
  }

  return 0;
}

int main(int argc, char **argv)
{
  std::string inputFileName;
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Help message")
    ("version,V", "Version number")
    ("input-file,i", po::value<std::string>(&inputFileName), "Input file name");

  try
  {
    po::positional_options_description p;
    p.add("input-file", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("version"))
    {
      std::cout << "Version: " << CSO_ANALYSIS_VER << std::endl;
      return 0;
    }
    if (vm.count("help"))
    {
      std::cout << "Usage: " << argv[0] << " [options]\n" << desc << std::endl;
      return 1;
    }

    if (!vm.count("input-file"))
    {
      std::cerr << "Missing input file" << std::endl;
      return 1;
    }

    return analyse(inputFileName.c_str());
  }
  catch (std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  return 0;
}

