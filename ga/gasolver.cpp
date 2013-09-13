/*
 * gasolver.h
 *
 *  Created on: 05-apr-2009
 *      Author: M. El-Kebir
 */

#include "gasolver.h"
#include "genotypetable.h"

GaSolver::GaSolver(const Data* pData, int nGenerations, int nRepetitions, double migrationProbability, bool verbose, bool storeFitnesses)
	: _pData(pData)
	, _nGenerations(nGenerations)
	, _nRepetitions(nRepetitions)
	, _pEliteIndividual(NULL)
	, _pOverAllBestIndividual(NULL)
	, _migrationProbability(migrationProbability)
	, _verbose(verbose)
	, _storeFitnesses(storeFitnesses)
{
	for (int i = 0; i < 4; i++)
	{
		_pBestIndividual[i] = NULL;
	}

	if (storeFitnesses)
	{
		for (int i = 0; i < 4000; i++)
		{
			_fitnessHistory.push_back("");
		}
	}
}

GaSolver::~GaSolver()
{
	delete _pEliteIndividual;
}

void GaSolver::generateInitialPopulation()
{
	for (int i = 0; i < 4; i++)
	{
		_population[i].clear();

		for (int j = 0; j < 1000; j++)
		{
			GaIndividual individual;
			_population[i].push_back(individual);
		}
	}
}

const GaIndividual& GaSolver::selectParent(int environment, int tournamentSize)
{
	assert(0 <= environment && environment < 4);

	GaIndividual* pResIndividual = NULL;
	double resFitness = DBL_MAX;

	for (int i = 0; i < tournamentSize; i++)
	{
		if (randDouble() < _migrationProbability)
		{
			int migEnv = (environment + randInt(4)) % 4;
			
			GaIndividual* pSelectedIndividual = &_population[migEnv][randInt(1000)];
			double selectedFitness = pSelectedIndividual->getFitness();

			if (selectedFitness < resFitness)
			{
				pResIndividual = pSelectedIndividual;
				resFitness = selectedFitness;
			}	
		}
		else
		{
			GaIndividual* pSelectedIndividual = &_population[environment][randInt(1000)];
			double selectedFitness = pSelectedIndividual->getFitness();

			if (selectedFitness < resFitness)
			{
				pResIndividual = pSelectedIndividual;
				resFitness = selectedFitness;
			}	
		}
	}

	assert(pResIndividual);
	return *pResIndividual;
}

void GaSolver::determineBestIndividuals()
{
	double overAllBestFitness = DBL_MAX;
	for (int i = 0; i < 4; i++)
	{
		double bestFitness = DBL_MAX;
		for (int j = 0; j < 1000; j++)
		{
			double individualFitness = _population[i][j].getFitness();
			if (individualFitness < bestFitness)
			{
				bestFitness = individualFitness;
				_pBestIndividual[i] = &_population[i][j];
			}
		}

		if (bestFitness < overAllBestFitness)
		{
			overAllBestFitness = bestFitness;
			_pOverAllBestIndividual = _pBestIndividual[i];
		}
	}

	double eliteFitness = _pEliteIndividual ? _pEliteIndividual->getFitness() : DBL_MAX;
	if (overAllBestFitness < eliteFitness)
	{
		delete _pEliteIndividual;
		_pEliteIndividual = new GaIndividual(*_pOverAllBestIndividual);
	}
}

void GaSolver::getTournamentSizes(const int environment, int& t1, int& t2)
{
	assert(0 <= environment && environment < 4);

	/** 
	 * environment 0: tournaments of 16 and 8
	 * environment 1: tournaments of 8 and 8
	 * environment 2: tournaments of 8 and 4
	 * environment 3: tournaments of 4 and 4
	 */

	switch (environment)
	{
	case 0:
		{
			t1 = 16;
			t2 = 8;
		}
		break;
	case 1:
		{
			t1 = 8;
			t2 = 8;
		}
		break;
	case 2:
		{
			t1 = 8;
			t2 = 4;
		}
		break;
	case 3:
		{
			t1 = 4;
			t2 = 4;
		}
		break;
	default:
		{
			// only needed to silence the compiler
			t1 = t2 = 16;
		}
		break;
	}
}

void GaSolver::solve()
{
	for (int l = 0; l < _nRepetitions; l++)
	{
		generateInitialPopulation();
		determineBestIndividuals();
		if (_storeFitnesses) updateFitnessHistory();

		for (int k = 0; k < _nGenerations; k++)
		{
			if (_verbose)
			{
				std::cout << "Iteration " << l << "\t" << "Generation " << k << "\t";
				for (int i = 0; i < 4; i++)
				{
					std::cout << _pBestIndividual[i]->getCost() << " " 
						<< _pBestIndividual[i]->getFitness() << "\t";
				}
				std::cout << GenotypeTable::getInstance()->getUseCount();
				std::cout << std::endl;
			}

			GaIndividualVector population[4];
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 250; j++)
				{
					int t1, t2;
					getTournamentSizes(i, t1, t2);
					GaIndividual parent1 = selectParent(i, t1);
					GaIndividual parent2 = selectParent(i, t2);

					population[i].push_back(parent1);
					population[i].push_back(parent2);

					population[i][4 * j].mutate();
					population[i][4 * j + 1].mutate();

					GaIndividual::crossover(parent1, parent2);

					population[i].push_back(parent1);
					population[i].push_back(parent2);
				}
			}

			for (int i = 0; i < 4; i++)
			{
				_population[i] = population[i];
			}

			determineBestIndividuals();
			if (_storeFitnesses) updateFitnessHistory();
		}
	}
}

void GaSolver::updateFitnessHistory()
{
	for (int i = 0; i < 1000; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			std::string& str = _fitnessHistory[i + j * 1000];
			if (str.compare(""))
			{
				str += ",";
			}
			
			double cost = _population[j][i].getCost();
			int badAllleles = _population[j][i].getNrOfBadAlleles();

			if (badAllleles == 0)
			{
				char buf[1024];
				sprintf(buf, "%f", cost);
				str += buf;
			}
			else
			{
				str += "\"\"";
			}
		}
	}
}

void GaSolver::printFitnessHistory(std::ostream& out) const
{
	for (int i = 0; i < 4000; i++)
	{
		out << _fitnessHistory[i].c_str() << std::endl;
	}
}
