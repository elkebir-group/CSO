/*
 * heuristicpurelinkage.h
 *
 *  Created on: 05-may-2009
 *      Author: M. El-Kebir
 */

#include "heuristicpurelinkage.h"
#include "../analysis/linkageanalysis.h"

const int HeuristicPureLinkage::_costMatrix[3][3] = { {0, 0, 0}, {10, 10, 1}, {10, 10, 10} };

HeuristicPureLinkage::HeuristicPureLinkage(int target, const LinkageAnalysis& linkageAnalysis)
	: HeuristicBase(target)
	, _linkageAnalysis(linkageAnalysis)
{
}

double HeuristicPureLinkage::getCost(const DpItem& item) const
{
	double cost = 0;
	for (int i = 0; i < _nLoci; i++)
	{
		for (int j = i + 1; j < _nLoci; j++)
		{
			LinkageType parentLinkage = _linkageAnalysis.getLinkageType(_target, i, j);
			// TODO: maak constistent ... (target eerste argument)
			LinkageType itemLinkage = item.getLinkage(i, j, _target);

			if (itemLinkage == Linked)
			{
				cost -= _linkageAnalysis.getLinkagePopPair(_target, i, j).second;
			}
			else if (itemLinkage == WeaklyLinked)
			{
				if (parentLinkage != Unlinked)
				{
					//cost *= RM[i][j];
					cost += _linkageAnalysis.getLinkagePopPair(_target, i, j).second;
				}
				else
				{
					//cost *= 0.5;
				}
			}
			else
			{
				//cost *= RM[i][j];
				cost += _linkageAnalysis.getLinkagePopPair(_target, i, j).second;
			}
			
			//cost += _costMatrix[itemLinkage][parentLinkage];
		}
	}

	return cost;
}
