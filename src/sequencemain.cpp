/*
 * sequencemain.cpp
 *
 *  Created on: 23-may-2018
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include <lemon/arg_parser.h>
#include <fstream>
#include "mutclonetree.h"

int main(int argc, char** argv)
{
  int seed = 0;
  int coverage;
  double purity = 1;
  int ploidy = 2;
  double seqErrorRate = 0;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("s", "Random number generator seed (default: 0)", seed, false)
    .refOption("C", "Coverage", coverage, true)
    .refOption("P", "Purity (default: 1)", purity, false)
    .refOption("ploidy", "Ploidy (default: 2)", ploidy, false)
    .refOption("E", "Per base sequencing error rate (default: 0)", seqErrorRate)
    .other("T", "Clone tree");
  ap.parse();
  
  MutCloneTree T;
  try
  {
    if (ap.files().size() >= 1)
    {
      std::string filenameTree = ap.files()[0];
      std::ifstream inT(filenameTree.c_str());
      if (!inT.good())
      {
        std::cerr << "Error: could not open '" << filenameTree << "' for reading" << std::endl;
        return 1;
      }
      inT >> T;
    }
    else
    {
      std::cin >> T;
    }
  }
  catch (std::runtime_error& e)
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  
  g_rng = std::mt19937(seed);
  std::cout << T.getReads(purity, coverage, seqErrorRate, ploidy);
}
