/*
 * preclustermain.cpp
 *
 *  Created on: 12-jun-2018
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "frequencymatrix.h"
#include <lemon/arg_parser.h>
#include <fstream>

int main(int argc, char** argv)
{
  std::string clusteringFilename;
  std::string frequencyFilename;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("C", "Clustering input filename", clusteringFilename)
    .other("F", "Frequency matrix");
  ap.parse();
  
  if (ap.files().size() != 1)
  {
    std::cerr << "Error: <F> must be specificed" << std::endl;
    return 1;
  }
  
  frequencyFilename = ap.files()[0];
  
  FrequencyMatrix F;
  try
  {
    if (frequencyFilename != "-")
    {
      std::ifstream inF(frequencyFilename.c_str());
      if (!inF.good())
      {
        std::cerr << "Error: failed to open '" << frequencyFilename << "' for reading" << std::endl;
        return 1;
      }
      
      inF >> F;
      inF.close();
    }
    else
    {
      std::cin >> F;
    }
  }
  catch (std::runtime_error& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  
  std::ifstream inC(clusteringFilename.c_str());
  if (!inC.good())
  {
    std::cerr << "Error: could not open '" << clusteringFilename << "' for reading" << std::endl;;
    return 1;
  }

  IntMatrix clustering = F.parseClustering(inC);
  inC.close();
  
  FrequencyMatrix newF = F.cluster(clustering);
  
  std::cout << newF;
  
  return 0;
}
