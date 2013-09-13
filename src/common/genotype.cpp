/*
 * genotype.cpp
 *
 *  Created on: 28-apr-2011
 *      Author: M. El-Kebir
 *
 */

#include "genotype.h"

std::istream& operator >>(std::istream& is, Genotype& genotype)
{
  return is >> genotype._c0 >> genotype._c1;
}

std::ostream& operator <<(std::ostream& os, const Genotype& genotype)
{
  return os << genotype._c0 << " " << genotype._c1;
}
