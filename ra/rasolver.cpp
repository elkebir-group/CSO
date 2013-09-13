/*
 * rasolver.cpp
 *
 *  Created on: 17-may-2009
 *      Author: M. El-Kebir
 */

#include "rasolver.h"
#include "coveranalysis.h"
#include "node.h"
#include "innernode.h"
#include "leafnode.h"

RaSolver::RaSolver(const int verbosity, const Data* pData, 
	const GenotypeTable* pGenotypeTable, int nRepetitions, int nTrialsPerCover)
	: _verbosity(verbosity)
	, _pData(pData)
	, _pGenotypeTable(pGenotypeTable)
	, _nRepetitions(nRepetitions)
	, _nTrialsPerCover(nTrialsPerCover)
{
	CoverAnalysis coverAnalysis(pData);
	GenotypeSet minimalCover = coverAnalysis.getMinimalCover();
	
	const int nLoci = _pData->getNumberOfLoci();
	_minimalCover = GenotypeList(minimalCover.begin(), minimalCover.end());
	_minimalCoverBitmap = std::vector<bool>(nLoci, false);

	const Genotype ideotype = _pData->getIdeotype();
	assert(ideotype.isHomozygous());

	for (GenotypeList::const_iterator it = _minimalCover.begin(); it != _minimalCover.end(); it++)
	{
		for (int i = 0; i < nLoci; i++)
		{
			// ideotype is homozygous
			_minimalCoverBitmap[i] = _minimalCoverBitmap[i] || 
				(*it)(0, i) == ideotype(0, i) || (*it)(1, i) == ideotype(0, i);
		}
	}
}

GenotypeList RaSolver::getParentCover()
{
	/* generate random cover */
	const int nLoci = _pData->getNumberOfLoci();
	GenotypeSet parents = _pData->getParents();
	const Genotype ideotype = _pData->getIdeotype();
	
	std::vector<bool> bitmap = _minimalCoverBitmap;

	bool done = true;
	/* check whether we are done */
	for (int i = 0; i < nLoci && done; i++)
	{
		if (!bitmap[i])
			done = false;
	}

	GenotypeList cover = _minimalCover;
	while (!done)
	{
		done = true;

		GenotypeSet::iterator it = parents.begin();
		int n = randInt((int) parents.size());
		for (int i = 0; i < n; i++) it++;
		cover.push_back(*it);

		for (int i = 0; i < nLoci; i++)
		{
			bitmap[i] = bitmap[i] || (*it)(0, i) == ideotype(0, i) || (*it)(1, i) == ideotype(0, i);;
			
			done &= bitmap[i];
		}

		parents.erase(it);
	}

	return cover;
}

GaIndividualPtr RaSolver::generateIndividual(const GenotypeList& cover)
{
	const Genotype& ideotype = _pData->getIdeotype();
	const int target = ideotype.getC0();
	const int nLoci = _pData->getNumberOfLoci();

	GenotypeNodeMap genotypeNodeMap;
	for (GenotypeList::const_iterator it = cover.begin(); it != cover.end(); it++)
	{
		const GenotypeGamete& parent = _pGenotypeTable->getGenotype(*it);
		
		/**
		 * Either parent._c0 or parent._c1 is propagated
		 * 
		 * Alternative: use parent.getGameteProportionate()._c; 
		 */
		int childGamete = randDouble() < 0.5 ? parent.getC0() : parent.getC1();
		NodePtr res = LeafNode::create(parent, false, childGamete, 1.0);
		
		genotypeNodeMap[parent] = res;
	}

	NodePtr homozygousGenotype;
	int cost = INT_MAX;
	NodePtr res = createCross(cover, genotypeNodeMap, target, homozygousGenotype, cost);
	if (homozygousGenotype && res->getGamete() != target)
	{
		/* backcross with homozygousGenotype if ideotype has not been reached yet */
		const GenotypeGamete& genotypeGamete = 
			_pGenotypeTable->getGenotype(homozygousGenotype->getGamete(), res->getGamete());
		
		res = InnerNode::create(res, homozygousGenotype, false, NoBackcross, 
			genotypeGamete.getGameteMostAlike(nLoci, target)._c);
	}
	if (res->getGenotype() != ideotype)
	{
		/* still no ideotype, then self */
		res = InnerNode::create(res, res, false, NoBackcross, res->getGamete());
	}

	return GaIndividualPtr(new GaIndividual(res));
}

NodePtr RaSolver::createCross(const GenotypeList& parents, GenotypeNodeMap& genotypeNodeMap, 
	int targetGamete, NodePtr& homozygousGenotype, int& costPopHomozygousGenotype)
{
	static const int nLoci = _pData->getNumberOfLoci();
	static const double gamma = _pData->getGamma();
	static const DoubleMatrix& RM = _pData->getRM();
	static const Genotype& ideotype = _pData->getIdeotype();

	int numberOfParents = (int) parents.size();
	if (numberOfParents == 1)
	{
		Genotype parent = parents.front();
		assert(genotypeNodeMap.find(parent) != genotypeNodeMap.end());
		
		NodePtr res = genotypeNodeMap[parent];
		
		int childGamete = res->getGamete();
		Genotype child = Genotype(childGamete, childGamete);
		
		/* self with probability 0.5, if not homozygous and offspring not already present */
		if (!parent.isHomozygous() && genotypeNodeMap.find(child) == genotypeNodeMap.end() && randDouble() < 0.5)
		{
			res = InnerNode::create(res, res, false, NoBackcross, res->getGamete());
			genotypeNodeMap[child] = res;
		}

		return res;
	}

	/* Partition parents into parents1 and parents2 at random */
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

	/* recurse on parents1 and parents2 */
	NodePtr p1 = createCross(parents1, genotypeNodeMap, targetGamete, homozygousGenotype, costPopHomozygousGenotype);
	NodePtr p2 = createCross(parents2, genotypeNodeMap, targetGamete, homozygousGenotype, costPopHomozygousGenotype);

	const Genotype& p1Genotype = p1->getGenotype();
	const Genotype& p2Genotype = p2->getGenotype();
	int p1Gamete = p1->getGamete();
	int p2Gamete = p2->getGamete();

	/* update homozygousGenotype if needed */
	if (p1Genotype.isHomozygous())
	{
		Genotype genotype(targetGamete, p1Gamete);
		int popCost = ideotype.computePop(nLoci, RM, gamma, genotype, genotype);

		/* check whether selfing (target, p1Gamete) is cheaper than costPopHomozygousGenotype */
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

		/* check whether selfing (target, p2Gamete) is cheaper than costPopHomozygousGenotype */
		if (popCost < costPopHomozygousGenotype)
		{
			costPopHomozygousGenotype = popCost;
			homozygousGenotype = p2;
		}
	}

	/* create new result node by crossing p1 and p2 */
	const GenotypeGamete& genotypeChild = _pGenotypeTable->getGenotype(p1Gamete, p2Gamete);

	NodePtr res;
	if (genotypeNodeMap.find(genotypeChild) == genotypeNodeMap.end())
	{
		const int gameteChild = genotypeChild.getGameteMostAlike(nLoci, targetGamete)._c;
		res = InnerNode::create(p1, p2, false, NoBackcross, gameteChild);
		genotypeNodeMap[genotypeChild] = res;
	}
	else
	{
		res = genotypeNodeMap[genotypeChild];
	}

	/* backcross with homozygousGenotype with probability 0.5 */
	if (homozygousGenotype && !res->getGenotype().isHomozygous() && randDouble() < 0.5)
	{
		const GenotypeGamete& genotypeGamete = 
			_pGenotypeTable->getGenotype(homozygousGenotype->getGamete(), res->getGamete());
		
		if (genotypeNodeMap.find(genotypeGamete) == genotypeNodeMap.end())
		{
			res = InnerNode::create(res, homozygousGenotype, false, NoBackcross, 
				genotypeGamete.getGameteMostAlike(nLoci, targetGamete)._c);
			genotypeNodeMap[genotypeGamete] = res;
		}
	}
	
	/* self with probability 0.2, if res is not homozygous */
	if (!res->getGenotype().isHomozygous() && randDouble() < 0.2)
	{
		Genotype selfedGenotype = Genotype(res->getGamete(), res->getGamete());
		if (genotypeNodeMap.find(selfedGenotype) == genotypeNodeMap.end())
		{
			res = InnerNode::create(res, res, false, NoBackcross, res->getGamete());
			genotypeNodeMap[selfedGenotype] = res;
		}
	}
	return res;
}

int RaSolver::solve()
{
	const Genotype& ideotype = _pData->getIdeotype();
	double bestCost = DBL_MAX;

	int correctIndividuals = 0;

	for (int i = 0; i < _nRepetitions; i++)
	{
		GenotypeList cover = getParentCover();

		for (int j = 0; j < _nTrialsPerCover; j++)
		{
			GaIndividualPtr pIndividual = generateIndividual(cover);
			if (pIndividual->getGenotype() == ideotype)
			{
				correctIndividuals++;

				double cost = pIndividual->getCost();

				if (_verbosity > 1)
				{
					printf("// Repetition %d. Trial %d. Generated correct individual with cost %f\n", i, j, cost);
				}

				if (cost < bestCost)
				{
					bestCost = cost;
					_pBestIndividual = pIndividual;

					if (_verbosity == 1)
					{
						printf("// Repetition %d. Trial %d. Generated new best individual with cost %f\n", i, j, cost);
					}
				}
			}
		}
	}

	return correctIndividuals;
}

std::string RaSolver::getMethodName() const
{
	char buf[1024];

	if (_pGenotypeTable->getRestrictGametes())
	{
		sprintf(buf, "RA-%d-%dg", _nRepetitions, _nTrialsPerCover);
	}
	else
	{
		sprintf(buf, "RA-%d-%d", _nRepetitions, _nTrialsPerCover);
	}

	return std::string(buf);
}
