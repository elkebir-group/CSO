/*
 * recipemain.cpp
 *
 *  Created on: 08-apr-2009
 *      Author: M. El-Kebir
 */

#define CSO_RECIPE_VERSION "08042009"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/program_options.hpp>
#include <exception>
#include "csodp.h"
#include "data.h"
#include "solver.h"
#include "dptablehashmap.h"
#include "dptablearray.h"
#include "dpitem.h"

namespace po = boost::program_options;

int recipeToDAG(const char* pInputFileName, const char* pRecipeFileName)
{
	Data* pData = Data::create(pInputFileName, false);
	if (!pData)	return 1;

	int nLoci = pData->getNumberOfLoci();
	const DoubleMatrix& RM = pData->getRM();
	double gamma = pData->getGamma();

	DpTableHashMap table(pData->getNumberOfLoci());
	Solver solver(pData, &table, 0, false);

	std::ifstream in(pRecipeFileName);
	if (!in.good())
	{
		std::cerr << "Failed to open file " << pRecipeFileName << std::endl;
		return 1;
	}

	char buf[1024];
	char* chromosomes[6];

	DpItem* pLastItem = NULL;
	int line = 1;
	while (!in.eof())
	{
		in.getline(buf, 1024);
		
		int i = 0;
		int j = 0;
		int n = strlen(buf);
		while (i++ < n && j < 6)
		{
			if (buf[i] == ' ' || buf[i] == 'x')
			{
				buf[i] = '\0';
				chromosomes[j++] = &buf[i+1];
			}
		}
		
		if (j == 0)
		{
			continue;
		}

		if (j != 6)
		{
			std::cerr << "Error in line " << line << std::endl;
			return 1;
		}

		int d0 = fromBitstring(nLoci, chromosomes[0]);
		int d1 = fromBitstring(nLoci, chromosomes[1]);
		int e0 = fromBitstring(nLoci, chromosomes[2]);
		int e1 = fromBitstring(nLoci, chromosomes[3]);
		int c0 = fromBitstring(nLoci, chromosomes[4]);
		int c1 = fromBitstring(nLoci, chromosomes[5]);

		Genotype C(c0, c1);
		Genotype D(d0, d1);
		Genotype E(e0, e1);

		DpItem* pItemD;
		DpItem* pItemE;
		
		if (table.isPresent(D)) 
			pItemD = &table.getItem(D);
		else 
			pItemD = &table.updateItem(DpItem(D, false));

		if (table.isPresent(E)) 
			pItemE = &table.getItem(E);
		else 
			pItemE = &table.updateItem(DpItem(E, false));
		

		if (table.isPresent(C))
		{
			std::cerr << "Error in line " << line << ". Genotype (";
			E.printGenotype(nLoci, false, std::cerr, "x");
			std::cerr << ") was already present" << std::endl;
			return 1;
		}

		unsigned long pop = C.computePop(nLoci, RM, gamma, D, E);
		pLastItem = &table.updateItem(DpItem(C, pItemD, pItemE, pop, 0, 0, 0, GenotypeSet(), false));

		line++;
	}
	
	pLastItem->updateAttributesFast(&table, true);
	solver.printDAG(*pLastItem, std::cout);

	delete pData;
	return 0;
}

int main(int argc, char **argv)
{
	std::string recipeFileName, inputFileName;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "Help message")
		("version,V", "Version number")
		("recipe-file,r", po::value<std::string>(&recipeFileName), "Recipe file name")
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
			std::cout << "Version: " << CSO_RECIPE_VERSION << std::endl;
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

		if (!vm.count("recipe-file"))
		{
			std::cerr << "Missing input file" << std::endl;
			return 1;
		}

		return recipeToDAG(inputFileName.c_str(), recipeFileName.c_str());
	}
	catch (...)
	{
		std::cerr << "Usage: " << argv[0] << " [options]\n" << desc << std::endl;
		return 1;
	}
}
