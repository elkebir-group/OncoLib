/*
 * mixmain.cpp
 *
 *  Created on: 22-may-2018
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include <lemon/arg_parser.h>
#include <fstream>
#include "mutclonetree.h"
#include "generatemixture.h"

int main(int argc, char** argv)
{
  int seed = 0;
  int nrSamplesPerAnatomicalSite = 0;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("s", "Random number generator seed (default: 0)", seed, false)
    .refOption("k", "Number of samples per anatomical site", nrSamplesPerAnatomicalSite, true)
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
  MutCloneTree sampledT = GenerateMixture(T).generate(nrSamplesPerAnatomicalSite);
  
  std::cout << sampledT;
  
  return 0;
}
