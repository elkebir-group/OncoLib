/*
 * tree2freqsmain.cpp
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
  lemon::ArgParser ap(argc, argv);
  ap.other("T", "Clone tree");
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
  
  std::cout << T.getFrequencies();
  
  return 0;
}
