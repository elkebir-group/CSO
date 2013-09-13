/* 
 * computeLB.cpp
 *
 *  Created on: 3-may-2011
 *      Author: M. El-Kebir
 */

#include <iostream>
#include <fstream>
#include <string>
#include <lemon/arg_parser.h>
#include "common/data.h"
#include "analysis/lowerbounds.h"

int main(int argc, char** argv)
{
  bool verbose = false;
  int nCrs;

  lemon::ArgParser ap(argc, argv);
  ap.refOption("v", "Verbose output", verbose, false)
    .refOption("c", "Number of crossings", nCrs, true);
  ap.parse();

  if (ap.files().size() == 0)
  {
    std::cerr << "Missing input file" << std::endl;
    return 1;
  }

  const std::string& inputFileName = ap.files()[0];
  Data* pData = Data::create(inputFileName.c_str());
  if (!pData) return 1;
  pData->updateGamma(nCrs);

  LowerBound LB(pData, verbose);
  std::cout << "Lower bounds" << std::endl
            << "============" << std::endl
            << "popMax: " << LB.getPopMaxLB() << std::endl
            << "pop: " << LB.getPopLB() << std::endl
            << "crs: " << LB.getCrossLB() << std::endl
            << "gen: " << LB.getGenLB() << std::endl;

  delete pData;
  return 0;
}
