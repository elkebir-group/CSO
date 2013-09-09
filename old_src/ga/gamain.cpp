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
#include "gasolver.h"
#include "stopwatch.h"

namespace po = boost::program_options;

int solve(int seed, int nGenerations, int nRepetitions, double migrationProbability, 
	const char* pInputFileName, bool useArray, bool verbose, bool allowPopMax, unsigned int limit,
	const char* pFitnessCsvFileName)
{
	Data* pData = Data::create(pInputFileName, true, allowPopMax);
	if (!pData)	return 1;

	const Genotype& ideotype = pData->getIdeotype();

	std::cout << "// Seed: " << seed << std::endl
		<< "// Migration probability: " << migrationProbability << std::endl
		<< "// Allow popMax: " << (allowPopMax ? "true" : "false") << std::endl;

	srand(seed);

	GenotypeTable* pTable = NULL;
	if (useArray)
	{
		if (verbose) std::cout << "\nAllocating array...\n";
		pTable = new GenotypeTableArray(limit);
	}
	else
	{
		pTable = new GenotypeTableHashMap(limit);
	}

	GaSolver solver(pData, nGenerations, nRepetitions, migrationProbability, verbose, strcmp(pFitnessCsvFileName, "") != 0);
	StopWatch stopWatch;
	stopWatch.start();
	solver.solve();
	stopWatch.stop();

	const GaIndividual* pEliteIndividual = solver.getEliteIndividual();
	if (pEliteIndividual->getGenotype() == ideotype)
	{
		pEliteIndividual->printDAG(std::cout);
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
			pInputFileName, "GA", nGenerations, nRepetitions, stopWatch.getElapsedTime(), 
			pEliteIndividual->getGen(), pEliteIndividual->getCross(), pEliteIndividual->getPop(),
			pEliteIndividual->getCost());
	}
	else
	{
		fprintf(stderr, "\"%s\",\"%s\",%d,%d,%.3f,\"-\",\"-\",\"-\",\"-\"\n",
			pInputFileName, "GA", nGenerations, nRepetitions, stopWatch.getElapsedTime());
	}

	if (strcmp(pFitnessCsvFileName, ""))
	{
		std::ofstream outFile(pFitnessCsvFileName);
		if (outFile.good())
		{
			solver.printFitnessHistory(outFile);
		}
		outFile.close();
	}

	delete pTable;
	delete pData;
	return 0;
}

int main(int argc, char **argv)
{
	unsigned int limit;
	bool useArray = false;
	bool verbose = false;
	bool allowPopMax = false;
	double migrationProbability;
	int nGenerations, nRepetitions;
	unsigned int seed;
	std::string inputFileName, fitnessCsvFileName;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "Help message")
		("version,V", "Version number")
		("verbose,v", po::value<bool>(&verbose)->zero_tokens(), "Verbose")
		("allow-pop-max", po::value<bool>(&allowPopMax)->zero_tokens(), "Allow popMax")
		("use-array,A", po::value<bool>(&useArray)->zero_tokens(), "Use array")
		("repetitions,r", po::value<int>(&nRepetitions)->default_value(10), "Number of repetitions")
		("generations,g", po::value<int>(&nGenerations)->default_value(100), "Number of generations per repetition")
		("migration-probability,m", po::value<double>(&migrationProbability)->default_value(0.001), "Migration probability")
		("seed,s" , po::value<unsigned int>(&seed)->default_value((unsigned int) time(NULL)), "Seed")
		("cache-size" , po::value<unsigned int>(&limit)->default_value(10000), "Cache size (#genotypes)")
		("fitness-csv", po::value<std::string>(&fitnessCsvFileName), "Fitness output file name (CSV)")
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

		solve(seed, nGenerations, nRepetitions, migrationProbability, 
			inputFileName.c_str(), useArray, verbose, allowPopMax, limit, fitnessCsvFileName.c_str());

		return 0;
	}
	catch (std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}
}
