/* 
 * processDAG.cpp
 *
 *  Created on: 28-apr-2011
 *      Author: M. El-Kebir
 */

#include <iostream>
#include <fstream>
#include <string>
#include <lemon/arg_parser.h>
#include "common/data.h"
#include "crossingschedule.h"

int main(int argc, char** argv)
{
  bool recomputePop;
  std::string dataFileName;

  lemon::ArgParser ap(argc, argv);
  ap.refOption("r", "Recompute population sizes", recomputePop, false)
    .refOption("d", "Data file name", dataFileName, true);
  ap.parse();

  if (ap.files().size() == 0)
  {
    std::cerr << "Missing input file" << std::endl;
    return 1;
  }

  const std::string& inputFileName = ap.files()[0];
  std::ifstream is(inputFileName.c_str());
  if (!is.good())
  {
    std::cerr << "Could not open input file" << std::endl;
    return 1;
  }

  Data* pData = Data::create(dataFileName.c_str());
  if (!pData) return 1;

  CrossingSchedule G(pData);
  G.loadDAG(is);

  pData->updateGamma(G.getCross());

  fprintf(stderr, "\"%s\",%lu,%lu,%.2f,%.2f,",
    dataFileName.c_str(),
    G.getGen(), G.getCross(), G.getPop(), G.getCost());

  if (recomputePop)
    G.recomputePop();
  G.printDAG(std::cout);

  fprintf(stderr, "%.3f,%.3f\n", G.getPop(), G.getCost());

  return 0;
}
