/*
 * gamain.cpp
 *
 *  Created on: 24-mar-2009
 *      Author: M. El-Kebir
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <boost/program_options.hpp>
#include <exception>
#include "csoga.h"
#include "data.h"
#include "genotypetablearray.h"
#include "genotypetablehashmap.h"
#include "gaindividual.h"
#include "leafnode.h"
#include "innernode.h"
#include "stopwatch.h"
#include "coveranalysis.h"
#include "rasolver.h"

namespace po = boost::program_options;

NodePtr createCross(const GenotypeList& parents, int targetGamete, NodePtr& homozygousGenotype, int& costPopHomozygousGenotype)
{
	static const int nLoci = Data::getInstance()->getNumberOfLoci();
	static const double gamma = Data::getInstance()->getGamma();
	static const DoubleMatrix& RM = Data::getInstance()->getRM();
	static const Genotype& ideotype = Data::getInstance()->getIdeotype();

	int numberOfParents = (int) parents.size();
	if (numberOfParents == 1)
	{
		const GenotypeGamete& parent = GenotypeTable::getInstance()->getGenotype(*parents.begin());
		int childGamete = parent.getGameteProportionate()._c;
		NodePtr res = LeafNode::create(parent, false, childGamete, 1.0);
		if (!res->getGenotype().isHomozygous() && randDouble() < 0.5)
		{
			res = InnerNode::create(res, res, false, NoBackcross, res->getGamete());
		}
		return res;
	}

	GenotypeList parents1, parents2 = parents;

	int parents1Size = 1 + randInt(numberOfParents - 1);
	for (int i = 0; i < parents1Size; i++)
	{
		GenotypeList::iterator it = parents2.begin();
		
		int n = randInt(parents2.size());
		for (int j = 0; j < n; j++) it++;
		
		parents1.push_back(*it);
		parents2.erase(it);
	}

	NodePtr p1 = createCross(parents1, targetGamete, homozygousGenotype, costPopHomozygousGenotype);
	NodePtr p2 = createCross(parents2, targetGamete, homozygousGenotype, costPopHomozygousGenotype);
	const Genotype& p1Genotype = p1->getGenotype();
	const Genotype& p2Genotype = p2->getGenotype();
	int p1Gamete = p1->getGamete();
	int p2Gamete = p2->getGamete();

	if (p1Genotype.isHomozygous())
	{
		Genotype genotype(targetGamete, p1Gamete);
		int popCost = ideotype.computePop(nLoci, RM, gamma, genotype, genotype);

		if (popCost < costPopHomozygousGenotype)
		{
			costPopHomozygousGenotype = popCost;
			homozygousGenotype = p1;
		}		
	}
	if (p2Genotype.isHomozygous())
	{
		Genotype genotype(targetGamete, p2Gamete);
		int popCost = ideotype.computePop(nLoci, RM, gamma, genotype, genotype);

		if (popCost < costPopHomozygousGenotype)
		{
			costPopHomozygousGenotype = popCost;
			homozygousGenotype = p2;
		}
	}

	const GenotypeGamete& genotypeChild = GenotypeTable::getInstance()->getGenotype(p1Gamete, p2Gamete);
	const int gameteChild = genotypeChild.getGameteMostAlike(nLoci, targetGamete)._c;

	NodePtr res = InnerNode::create(p1, p2, false, NoBackcross, gameteChild);
	if (homozygousGenotype && !res->getGenotype().isHomozygous() && randDouble() < 0.5)
	{
		const GenotypeGamete& genotypeGamete = 
			GenotypeTable::getInstance()->getGenotype(homozygousGenotype->getGamete(), gameteChild);
		
		res = InnerNode::create(res, homozygousGenotype, false, NoBackcross, 
			genotypeGamete.getGameteMostAlike(nLoci, targetGamete)._c);
	}
	if (!res->getGenotype().isHomozygous() && randDouble() < 0.5)
	{
		res = InnerNode::create(res, res, false, NoBackcross, res->getGamete());
	}
	return res;
}

GaIndividualPtr generateIndividual(const Data* pData)
{
	GenotypeList parents(pData->getParents().begin(), pData->getParents().end());
	const Genotype& ideotype = Data::getInstance()->getIdeotype();

	/* generate random cover */
	int nLoci = pData->getNumberOfLoci();
	std::vector<bool> bitmap[2];
	bitmap[0] = bitmap[1] = std::vector<bool>(nLoci, false);

	GenotypeList cover;
	bool done = false;
	while (!done)
	{
		done = true;

		GenotypeList::iterator it = parents.begin();
		int n = randInt((int) parents.size());
		for (int i = 0; i < n; i++) it++;
		cover.push_back(*it);

		for (int i = 0; i < nLoci; i++)
		{
			bitmap[0][i] = bitmap[0][i] || (*it)(0, i) == ideotype(0, i);
			bitmap[1][i] = bitmap[1][i] || (*it)(1, i) == ideotype(1, i);
			
			done &= (bitmap[0][i] || bitmap[1][i]);
		}

		parents.erase(it);
	}

	NodePtr homozygousGenotype;
	int cost = INT_MAX;
	NodePtr res = createCross(cover, ideotype.getC0(), homozygousGenotype, cost);
	if (homozygousGenotype && res->getGamete() != ideotype.getC0())
	{
		const GenotypeGamete& genotypeGamete = 
			GenotypeTable::getInstance()->getGenotype(homozygousGenotype->getGamete(), res->getGamete());
		
		res = InnerNode::create(res, homozygousGenotype, false, NoBackcross, 
			genotypeGamete.getGameteMostAlike(nLoci, ideotype.getC0())._c);
	}
	/* TODO: add backcross */
	if (res->getGenotype() != ideotype)
	{
		res = InnerNode::create(res, res, false, NoBackcross, res->getGamete());
	}

	return GaIndividualPtr(new GaIndividual(res));
}

int solve(int seed, int nRepetitions, int nTrials, const char* pInputFileName, 
	bool useArray, int verbosity, unsigned int limit, int maxCrossOver, bool restrictGametes)
{
	Data* pData = Data::create(pInputFileName);
	if (!pData)	return 1;

	if (!pData->getIdeotype().isHomozygous())
	{
		std::cout << "This method only works for homozygous ideotypes" << std::endl;
		return 1;
	}

	std::cout << "// Seed: " << seed << std::endl; 
	srand(seed);

	GenotypeTable* pTable = NULL;
	if (useArray)
	{
		std::cout << "\nAllocating array...\n";
		pTable = new GenotypeTableArray(limit, maxCrossOver, restrictGametes);
	}
	else
	{
		pTable = new GenotypeTableHashMap(limit, maxCrossOver, restrictGametes);
	}

	StopWatch stopWatch;
	stopWatch.start();
	
	RaSolver raSolver(verbosity, pData, pTable, nRepetitions, nTrials);
	int correctIndividuals = raSolver.solve();
	GaIndividualPtr pBestIndividual = raSolver.getBestIndividual();

	stopWatch.stop();

	if (verbosity > 0)
	{
		std::cout << "// Number of considered genotypes: " << pTable->getUseCount() << std::endl;
		std::cout << "// Number of correct individuals: " << correctIndividuals << std::endl;
	}

	if (pBestIndividual)
	{	
		pBestIndividual->printDAG(std::cout);

		/* Format:
		 * - Input file name
		 * - Method name
		 * - Method specifics 1
		 * - Method specifics 2
		 * - Elapsed process time
		 * - Number of generations
		 * - Number of crossings
		 * - Population size
		 * - Total cost
		 */
		fprintf(stderr, "\"%s\",\"%s\",%d,%d,%.3f,%lu,%lu,%lu,%.2f\n",
			pInputFileName, raSolver.getMethodName().c_str(), nRepetitions, nTrials, stopWatch.getElapsedTime(), 
			pBestIndividual->getGen(), pBestIndividual->getCross(), pBestIndividual->getPop(),
			pBestIndividual->getCost());
	}
	else
	{
		fprintf(stderr, "\"%s\",\"%s\",%d,%d,%.3f,\"-\",\"-\",\"-\",\"-\"\n",
			pInputFileName, raSolver.getMethodName().c_str(), nRepetitions, nTrials, stopWatch.getElapsedTime());
	}

	delete pTable;
	delete pData;
	return pBestIndividual ? true : false;
}

int main(int argc, char **argv)
{
	unsigned int limit;
	bool useArray = false;
	int verbosity = 0;
	int nRepetitions;
	int nTrials;
	unsigned int seed;
	bool restrictGametes = false;
	int maxCrossOver;
	std::string inputFileName;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "Help message")
		("version,V", "Version number")
		("verbosity-level,v", po::value<int>(&verbosity)->default_value(0), "verbosity level")
		("use-array,A", po::value<bool>(&useArray)->zero_tokens(), "Use array")
		("repetitions,r", po::value<int>(&nRepetitions)->default_value(1000), "Number of repetitions")
		("trials,t", po::value<int>(&nTrials)->default_value(10), "Number of trials per repetition")
		("seed,s" , po::value<unsigned int>(&seed)->default_value((unsigned int) time(NULL)), "Seed")
		("cache-size" , po::value<unsigned int>(&limit)->default_value(10000), "Cache size (#genotypes)")
		("max-crossover,c", po::value<int>(&maxCrossOver)->default_value(-1), "Maximum number of crossovers")
		("restrict-gametes,g", po::value<bool>(&restrictGametes)->zero_tokens(), "Restrict gametes")
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
			std::cout << "Version: " << CSO_GA_VERSION << std::endl;
			return 0;
		}
		if (vm.count("help"))
		{
			std::cout << "Usage: " << argv[0] << " [options]\n" << desc << std::endl;
			return 0;
		}

		if (!vm.count("input-file"))
		{
			std::cerr << "Missing input file" << std::endl;
			return 1;
		}

		solve(seed, nRepetitions, nTrials, inputFileName.c_str(), useArray, verbosity, limit, maxCrossOver, restrictGametes);

		return 0;
	}
	catch (std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}
}
