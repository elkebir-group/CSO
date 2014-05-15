/*
 * heuristics.h
 *
 *  Created on: 03-apr-2009
 *      Author: s030858
 */

#include <functional>
#include "cso.h"
#include "data.h"
#include "dpitem.h"

struct Heuristic1Predicate : public std::binary_function<DpItem, DpItem, bool>
{
  bool operator()(const DpItem& item1, const DpItem& item2) const
  {
    const Genotype& ideotype = Data::getInstance()->getIdeotype();
    int nLoci = Data::getInstance()->getNumberOfLoci();

    int ideotype0 = ideotype.getC0();
    int ideotype1 = ideotype.getC1();
    int item1_0 = item1.getC0();
    int item1_1 = item1.getC1();
    int item2_0 = item2.getC0();
    int item2_1 = item2.getC1();

    int correctAlleles1 = 0, correctAlleles2 = 0;
    for (int i = 0; i < nLoci; i++)
    {
      ideotype0 >>= 1;
      ideotype1 >>= 1;
      item1_0 >>= 1;
      item1_1 >>= 1;
      item2_0 >>= 1;
      item2_1 >>= 1;

      if ((ideotype0 & 1) == (item1_0 & 1) || (ideotype0 & 1) == (item1_1 & 1) ||
        (ideotype1 & 1) == (item1_0 & 1) || (ideotype1 & 1) == (item1_1 & 1))
      {
        correctAlleles1++;
      }

      if ((ideotype0 & 1) == (item2_0 & 1) || (ideotype0 & 1) == (item2_1 & 1) ||
        (ideotype1 & 1) == (item2_0 & 1) || (ideotype1 & 1) == (item2_1 & 1))
      {
        correctAlleles2++;
      }
    }

    return correctAlleles1 > correctAlleles2;
  }
};

struct Heuristic2Predicate : public std::binary_function<DpItem, DpItem, bool>
{
  bool operator()(const DpItem& item1, const DpItem& item2) const
  {
    const Genotype& ideotype = Data::getInstance()->getIdeotype();
    int nLoci = Data::getInstance()->getNumberOfLoci();

    int ideotype0 = ideotype.getC0();
    int ideotype1 = ideotype.getC1();
    int item1_0 = item1.getC0();
    int item1_1 = item1.getC1();
    int item2_0 = item2.getC0();
    int item2_1 = item2.getC1();

    int correctAlleles1 = std::max(nLoci - numberOfDifferences(nLoci, item1_0, ideotype0), nLoci - numberOfDifferences(nLoci, item1_1, ideotype0)) +
      std::max(nLoci - numberOfDifferences(nLoci, item1_0, ideotype1), nLoci - numberOfDifferences(nLoci, item1_1, ideotype1));
    int correctAlleles2 = std::max(nLoci - numberOfDifferences(nLoci, item2_0, ideotype0), nLoci - numberOfDifferences(nLoci, item2_1, ideotype0)) +
      std::max(nLoci - numberOfDifferences(nLoci, item2_0, ideotype1), nLoci - numberOfDifferences(nLoci, item2_1, ideotype1));

    return correctAlleles1 > correctAlleles2;
  }
};

struct Heuristic3Predicate : public std::binary_function<DpItem, DpItem, bool>
{
  bool operator()(const DpItem& item1, const DpItem& item2) const
  {
    double cost1 = Data::getInstance()->getCost(item1.getCumPop(), item1.getGen(), item1.getCumCross());
    double cost2 = Data::getInstance()->getCost(item2.getCumPop(), item2.getGen(), item2.getCumCross());
    return cost1 < cost2;
  }
};

struct Heuristic4Predicate : public std::binary_function<DpItem, DpItem, bool>
{
private:
  int compute(int nLoci, const DoubleMatrix& RM, double gamma, const DpItem& item, int target) const
  {
    //std::cout << "//";
    //item.printGenotype(nLoci, false);
    //std::cout << "\t" << target;
    int c0 = item.getC0();
    int c1 = item.getC1();
    int penalty = 1;

    int prevIndex = -1;
    double prob = 1;
    for (int i = 0; i < nLoci; i++)
    {
            int c0_bit = GET_BIT(nLoci, c0, i);
            int c1_bit = GET_BIT(nLoci, c1, i);
            int target_bit = GET_BIT(nLoci, target, i);

      if (c0_bit != target_bit && c1_bit != target_bit)
      {
        penalty++;
        // flip target_bit
        if (target_bit == 1)
                    target &= ~(1 << (nLoci - i - 1));
        else
                    target |= (1 << (nLoci - i - 1));

        if (prevIndex != -1)
        {
          prob *= (1 - RM[prevIndex][i]);
          prevIndex = i;
        }
        else
        {
          prob = 0.5;
          prevIndex = i;
        }
      }
    }
    //std::cout << "\t" << target << "\t" << probToPop(item.computeProb(nLoci, RM, target), gamma) << "\t" << penalty * 100  << std::endl;
    return probToPop(item.computeProb(nLoci, RM, target), gamma) + penalty * 100;// + probToPop(prob, gamma);
  }

public:
  bool operator()(const DpItem& item1, const DpItem& item2) const
  {
    int nLoci = Data::getInstance()->getNumberOfLoci();
    const DoubleMatrix& RM = Data::getInstance()->getRM();
    double gamma = Data::getInstance()->getGamma();
    const Genotype& ideotype = Data::getInstance()->getIdeotype();

    double cost1 = 0;//Data::getInstance()->getCost(item1.getCumPop(), item1.getGen(), item1.getCumCross());
    double cost2 = 0;//Data::getInstance()->getCost(item2.getCumPop(), item2.getGen(), item2.getCumCross());

    cost1 += std::min(compute(nLoci, RM, gamma, item1, ideotype.getC0()), compute(nLoci, RM, gamma, item1, ideotype.getC1()));
    cost2 += std::min(compute(nLoci, RM, gamma, item2, ideotype.getC0()), compute(nLoci, RM, gamma, item2, ideotype.getC1()));

    return cost1 < cost2;
  }
};
