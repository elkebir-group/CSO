/*
 * dpmain.cpp
 *
 *  Created on: 17-feb-2009
 *      Author: M. El-Kebir
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/program_options.hpp>
#include <exception>
#include "csodp.h"
#include "data.h"
#include "solver.h"
#include "heuristicsolver.h"
#include "dptablehashmap.h"
#include "dptablearray.h"
#include "dpitem.h"
#include "stopwatch.h"

namespace po = boost::program_options;

int solve(int maxIt, const char* pInputFileName, bool useArray, 
	int verbosity, int heuristic, int maxNewItemSetSize, 
	int linkageAnalysis, int maxCrossOver, 
	bool restrictGametes, bool restrictSelfing)
{
	Data* pData = Data::create(pInputFileName);
	if (!pData)	return 1;

    pData->printRM();

	const Genotype& ideotype = pData->getIdeotype();

	DpTable* pTable = NULL;
	if (useArray)
	{
		if (verbosity > 0) 
			std::cout << "\nAllocating array...\n";
		pTable = new DpTableArray(maxCrossOver, restrictGametes);
	}
	else
	{
		pTable = new DpTableHashMap(maxCrossOver, restrictGametes);
	}

	Solver* pDP = NULL;
	if (heuristic)
	{
		pDP = new HeuristicSolver(pData, pTable, maxIt, verbosity, (HeuristicType) heuristic, maxNewItemSetSize, linkageAnalysis, restrictSelfing);
	}
	else
	{
		pDP = new Solver(pData, pTable, maxIt, verbosity);
	}

	StopWatch stopWatch;
	stopWatch.start();
	pDP->solve();
	stopWatch.stop();
	
	if (pTable->isPresent(ideotype))
	{
		//std::cout << pDP->getResDAG();

		const DpItem& ideotypeItem = pTable->getItem(ideotype);
		pDP->printDAG(ideotypeItem, std::cout);
		
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
		if (heuristic != 0)
		{
			fprintf(stderr, "\"%s\",\"%s\",%d,%d,%.3f,%lu,%lu,%lu,%.2f\n",
				pInputFileName, pDP->getMethodName().c_str(), heuristic, maxNewItemSetSize, stopWatch.getElapsedTime(), 
				ideotypeItem.getGen(), ideotypeItem.getCumCross(), ideotypeItem.getCumPop(),
				pDP->getScore(ideotypeItem));
		}
		else
		{
			fprintf(stderr, "\"%s\",\"%s\",\"-\",\"-\",%.3f,%lu,%lu,%lu,%.2f\n",
				pInputFileName, pDP->getMethodName().c_str(), stopWatch.getElapsedTime(), 
				ideotypeItem.getGen(), ideotypeItem.getCumCross(), ideotypeItem.getCumPop(),
				pDP->getScore(ideotypeItem));
		}
	}
	else
	{
		if (heuristic != 0)
		{
			fprintf(stderr, "\"%s\",\"%s\",%d,%d,%.3f,\"-\",\"-\",\"-\",\"-\"\n",
				pInputFileName, pDP->getMethodName().c_str(), heuristic, maxNewItemSetSize, stopWatch.getElapsedTime());
		}
		else
		{
			fprintf(stderr, "\"%s\",\"%s\",\"-\",\"-\",%.3f,\"-\",\"-\",\"-\",\"-\"\n",
				pInputFileName, pDP->getMethodName().c_str(), stopWatch.getElapsedTime());
		}
	}
	
	delete pTable;
	delete pData;
	
	return 0;
}

int main(int argc, char **argv)
{
	int linkageAnalysis;
	bool useArray = false;
	int verbosity = 0;
	int maxIt;
	int heuristic, maxNewItemSize, maxCrossOver;
	bool restrictGametes = false, restrictSelfing = false;
	std::string inputFileName;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "Help message")
		("version,V", "Version number")
		("verbosity-level,v", po::value<int>(&verbosity)->default_value(0), "verbosity level")
		("use-array,A", po::value<bool>(&useArray)->zero_tokens(), "Use array")
		("max-iterations,m", po::value<int>(&maxIt)->default_value(INT_MAX), "Maximum number of iterations")
		("max-genotypes,n", po::value<int>(&maxNewItemSize)->default_value(50), "Maximum number of new genotypes per iteration (requires heuristics)")
		("linkage-analysis,l", po::value<int>(&linkageAnalysis)->default_value(0), "Maximum number of new genotypes due to linkage analysis (requires heuristics)")
		("heuristic,h", po::value<int>(&heuristic)->default_value(0), "Heuristic option")
		("max-crossover,c", po::value<int>(&maxCrossOver)->default_value(-1), "Maximum number of crossovers")
		("restrict-gametes,g", po::value<bool>(&restrictGametes)->zero_tokens(), "Restrict gametes")
		("restrict-selfing,s", po::value<bool>(&restrictSelfing)->zero_tokens(), "Restrict selfing (no intermediary selfing, expect for P)")
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
			std::cout << "Version: " << CSO_DP_VERSION << std::endl;
			return 0;
		}
		if (vm.count("help"))
		{
			std::cout << "Usage: " << argv[0] << " [options]\n" << desc << std::endl;
			return 0;
		}

		if (vm.count("heuristic"))
		{
			if (!(0 <= heuristic && heuristic < 13))
			{
				std::cerr << "Value of --heuristic must be between 1 and 12" << std::endl;
				return 1;
			}
		}

		if (!vm.count("input-file"))
		{
			std::cerr << "Missing input file" << std::endl;
			return 1;
		}

		return solve(maxIt, inputFileName.c_str(), useArray, verbosity, 
			heuristic, maxNewItemSize, linkageAnalysis, 
			maxCrossOver, restrictGametes, restrictSelfing);
	}
	catch (...)
	{
		std::cerr << "Usage: " << argv[0] << " [options]\n" << desc << std::endl;
		return 1;
	}
}
