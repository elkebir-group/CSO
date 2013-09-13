/*
 * gasolver.h
 *
 *  Created on: 05-apr-2009
 *      Author: M. El-Kebir
 */

#ifndef GASOLVER_H_
#define GASOLVER_H_

#include "cso.h"
#include "csoga.h"
#include "data.h"
#include "gaindividual.h"

class GaSolver
{
private:
	const Data* _pData;
	const int _nGenerations;
	const int _nRepetitions;
	GaIndividualVector _population[4];
	/* Most elite individual among all generations and repetitions */
	GaIndividual* _pEliteIndividual;
	/* Overall best individual among all environments */
	GaIndividual* _pOverAllBestIndividual;
	/* Best individual per environment */
	GaIndividual* _pBestIndividual[4];
	const double _migrationProbability;
	const bool _verbose;
	const bool _storeFitnesses;
	std::vector<std::string> _fitnessHistory;
	
	void generateInitialPopulation();
	const GaIndividual& selectParent(int environment, int tournamentSize);
	void getTournamentSizes(const int environment, int& t1, int& t2);
	void determineBestIndividuals();
	void updateFitnessHistory();

public:
	GaSolver(const Data* pData, int nGenerations, int nRepetitions, double migrationProbability, bool verbose, bool storeFitnesses);
	~GaSolver();
	void solve();
	const GaIndividual* getEliteIndividual() const;
	void printFitnessHistory(std::ostream& out) const;
};

inline const GaIndividual* GaSolver::getEliteIndividual() const
{
	return _pEliteIndividual;
}

#endif
