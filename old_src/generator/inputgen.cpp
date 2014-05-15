/*
 * inputgen.cpp
 *
 *  Created on: 31-mar-2009
 *      Author: M. El-Kebir
 */

#define CSO_GEN_VER "01042009"

#include <stdio.h>
#include <iostream>
#include <boost/program_options.hpp>
#include <exception>
#include <time.h>
#include "genotype.h"
#include "cso.h"
#include "linkageanalysis.h"

namespace po = boost::program_options;

int getUncoveredLocusIndex(int nLoci, const GenotypeSet& parents, const Genotype& ideotype)
{
  unsigned int res = 0;
  for (GenotypeSet::const_iterator it = parents.begin(); it != parents.end(); it++)
  {
    res |= it->getC0() | it->getC1();
  }

  res = (~res) & ((1 << nLoci) - 1);
  for (int i = 0; i < nLoci; i++)
  {
    if (res & 1) return i;
        res >>= 1;
    }

  return nLoci;
}

void generateParents(GenotypeSet& parents, int nLoci, int nParents, int maxCorrectAlleles, double homozygousProb, const Genotype& ideotype)
{
  while ((int) parents.size() != nParents)
  {
    int c0 = randChromosome(nLoci);
    while (numberOfDifferences(nLoci, c0, ideotype.getC0()) < (nLoci - maxCorrectAlleles) ||
      numberOfDifferences(nLoci, c0, ideotype.getC1()) < (nLoci - maxCorrectAlleles))
    {
      c0 = randChromosome(nLoci);
    }
    int c1;

    if (randDouble() < homozygousProb)
    {
      c1 = c0;
    }
    else
    {
      c1 = randChromosome(nLoci);
      while (numberOfDifferences(nLoci, c1, ideotype.getC0()) < (nLoci - maxCorrectAlleles) ||
        numberOfDifferences(nLoci, c1, ideotype.getC1()) < (nLoci - maxCorrectAlleles))
      {
        c1 = randChromosome(nLoci);
      }
    }

    Genotype genotype(c0, c1);
    parents.insert(genotype);
  }

  /* make sure parents can lead to ideotype */
  int index;
  while ((index = getUncoveredLocusIndex(nLoci, parents, ideotype)) != nLoci)
  {
    int n = 1 + randInt(nParents / 10);
    for (int i = 0; i < n; i++)
    {
      int k = randInt(nParents);
      GenotypeSet::iterator it = parents.begin();
      for (int j = 0; j < k; j++) it++;

      int c0 = it->getC0(), c1 = it->getC1();
      if (c0 == c1)
      {
        c0 |= (1 << index);
        c1 = c0;
      }
      else if (randDouble() < 0.5)
      {
        c0 |= (1 << index);
      }
      else
      {
        c1 |= (1 << index);
      }

      Genotype alteredParent(c0, c1);
      parents.erase(it);
      parents.insert(alteredParent);
    }
  }
}
void generate(unsigned int seed, int nLoci, int nParents, int maxCorrectAlleles, double gamma,
    int costGen, int costCrossOver, int costNode, int popMax,
  const Genotype& ideotype, const DoubleVector& mapDistances, double homozygousProb)
{
  char buf1[1024];
  char buf2[1024];

  const DoubleMatrix RM = Data::generateRM(mapDistances);

  GenotypeSet parents;
  generateParents(parents, nLoci, nParents, maxCorrectAlleles, homozygousProb, ideotype);

  LinkageAnalysis analysis(nLoci, parents, ideotype, RM, gamma);
  unsigned long inferredPopMax = analysis.getPopMax();
  if ((unsigned long) popMax < inferredPopMax)
    popMax = (int) (1.5 * inferredPopMax);

  printf("<!-- seed: %u, max number of correct alleles: %d -->\n", seed, maxCorrectAlleles);
    printf("<CSO nLoci=\"%d\" gamma=\"%f\" popMax=\"%d\" costCrossOver=\"%d\" costGen=\"%d\" costNode=\"%d\">\n",
        nLoci, gamma, popMax, costCrossOver, costGen, costNode);

  printf("\t<Parents>\n");
  for (GenotypeSet::const_iterator it = parents.begin(); it != parents.end(); it++)
  {
    toBitstring(it->getC0(), nLoci, buf1);
    toBitstring(it->getC1(), nLoci, buf2);

    printf("\t\t<Parent c0=\"%s\" c1=\"%s\"/>\n", buf1, buf2);
  }
  printf("\t</Parents>\n");

  toBitstring(ideotype.getC0(), nLoci, buf1);
  toBitstring(ideotype.getC1(), nLoci, buf2);
  printf("\t<Ideotype c0=\"%s\" c1=\"%s\"/>\n", buf1, buf2);
  printf("\t<cM>");
  for (std::vector<double>::const_iterator it = mapDistances.begin(); it != mapDistances.end(); it++)
  {
    printf("%.2f ", *it);
  }
  printf("</cM>\n");
  printf("</CSO>\n");
}

int main(int argc, char **argv)
{
  int nLoci = -1, nParents = -1, maxCorrectAlleles;
  unsigned int seed;
  std::string ideotypeString;
  std::vector<double> mapDistances;
  double homozygousProb;
    int costGen, costCrossOver, costNode, popMax;
  double gamma;
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Help message")
    ("version,v", "Version number")
        ("seed,s" , po::value<unsigned int>(&seed)->default_value((unsigned int) time(NULL)), "Seed")
    ("number-loci,m", po::value<int>(&nLoci), "Number of loci")
    ("number-parents,n", po::value<int>(&nParents), "Number of parents")
    ("max-correct-alleles,a", po::value<int>(&maxCorrectAlleles), "Maximum number of correct alleles")
    ("homozygous-prob,p", po::value<double>(&homozygousProb)->default_value(0.5), "Probability of generating a homozygous parent")
    ("gamma", po::value<double>(&gamma)->default_value(0.99), "Gamma")
    ("pop-max", po::value<int>(&popMax)->default_value(500), "Maximum population size")
    ("cost-gen", po::value<int>(&costGen)->default_value(100), "Cost of a generation")
        ("cost-cross", po::value<int>(&costNode)->default_value(100), "Cost of a crossing")
        ("cost-pop", po::value<int>(&costCrossOver)->default_value(1), "Cost of an individual")
    ("ideotype,i", po::value<std::string>(&ideotypeString), "Ideotype")
    ("map-distances", po::value<std::vector<double> >(&mapDistances)->composing(), "Map distances in cM");

  try
  {
    po::positional_options_description p;
    p.add("map-distances", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    srand(seed);

    if (vm.count("version"))
    {
            std::cout << "Version: " << CSO_GEN_VER << std::endl;
      return 0;
    }
    if (vm.count("help"))
    {
      std::cout << "Usage: " << argv[0] << " [options]\n" << desc << std::endl;
      return 0;
    }

    if (!vm.count("number-loci"))
    {
      std::cerr << "Missing number of loci" << std::endl;
      return 1;
    }

    if (!vm.count("number-parents"))
    {
      std::cerr << "Missing number of parents" << std::endl;
      return 1;
    }

    if (!vm.count("max-correct-alleles"))
    {
      std::cerr << "Missing maximum number of correct alleles" << std::endl;
      return 1;
    }
    else if (maxCorrectAlleles > nLoci)
    {
      std::cerr << "Invalid value of maximum number of correct alleles" << std::endl;
      return 1;
    }

    Genotype ideotype((1 << nLoci) - 1, (1 << nLoci) - 1);
    if (vm.count("ideotype"))
    {
      if ((int) ideotypeString.size() != 2 * nLoci + 1)
      {
        std::cerr << "Invalid ideotype" << std::endl;
        return 1;
      }
      else
      {
        const char* pIdeotypeStr = ideotypeString.c_str();
        int c0 = fromBitstring(nLoci, pIdeotypeStr);
        int c1 = fromBitstring(nLoci, pIdeotypeStr + nLoci + 1);

        ideotype = Genotype(c0, c1);
      }
    }

    DoubleMatrix RM(nLoci);
    for (int i = 0; i < nLoci; i++)
    {
      RM[i] = DoubleVector(nLoci);
    }

    if (mapDistances.size())
    {
      if ((int) mapDistances.size() != nLoci)
      {
        std::cerr << "Expected " << nLoci << " map distances." << std::endl;
        return 1;
      }
    }
    else
    {
      mapDistances = std::vector<double>(nLoci);
      mapDistances[0] = randInt(15);
      for (int i = 1; i < nLoci; i++)
      {
        mapDistances[i] += mapDistances[i - 1] + 1.5 + (randInt(30) / 2.0);
      }
    }

        generate(seed, nLoci, nParents, maxCorrectAlleles, gamma, costGen, costCrossOver, costNode, popMax, ideotype, mapDistances, homozygousProb);

    return 0;
  }
  catch (std::exception& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}

