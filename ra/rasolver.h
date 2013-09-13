/*
 * rasolver.h
 *
 *  Created on: 17-may-2009
 *      Author: M. El-Kebir
 */

#ifndef RASOLVER_H_
#define RASOLVER_H_

#include "cso.h"
#include "csoga.h"
#include "data.h"
#include "genotypetable.h"
#include <map>

typedef std::map<Genotype, NodePtr> GenotypeNodeMap;

class RaSolver
{
private:
	const int _verbosity;
	const Data* _pData;
	const GenotypeTable* _pGenotypeTable;
	const int _nRepetitions;
	const int _nTrialsPerCover;
	GenotypeList _minimalCover;
	std::vector<bool> _minimalCoverBitmap;

	GaIndividualPtr _pBestIndividual;
	
	GenotypeList getParentCover();
	NodePtr createCross(const GenotypeList& parents, GenotypeNodeMap& genotypeNodeMap,
		int targetGamete, NodePtr& homozygousGenotype, int& costPopHomozygousGenotype);
	GaIndividualPtr generateIndividual(const GenotypeList& cover);

public:
	RaSolver(const int verbosity, const Data* pData, 
		const GenotypeTable* pGenotypeTable, int nRepetitions, int nTrialsPerCover);
	int solve();
	GaIndividualPtr getBestIndividual() const;
	std::string getMethodName() const;
};

inline GaIndividualPtr RaSolver::getBestIndividual() const
{
	return _pBestIndividual;
}

#endif /* RASOLVER_H_ */