/*
 * reads2freqsmain.cpp
 *
 *  Created on: 23-may-2018
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include <lemon/arg_parser.h>
#include <fstream>
#include "readmatrix.h"

int main(int argc, char** argv)
{
  double alpha = 0.01;
  lemon::ArgParser ap(argc, argv);
  ap.refOption("a", "Confidence interval (default: 0.01)", alpha)
    .other("R", "Read matrix");
  ap.parse();
  
  if (ap.files().size() != 1)
  {
    std::cerr << "Error: <R> must be specificed" << std::endl;
    return 1;
  }
  
  std::string filenameR = ap.files()[0].c_str();
  ReadMatrix R;
  try
  {
    if (filenameR != "-")
    {
      std::ifstream inR(filenameR.c_str());
      if (!inR.good())
      {
        std::cerr << "Error: failed to open '" << ap.files()[0] << "' for reading" << std::endl;
        return 1;
      }
      
      inR >> R;
    }
    else
    {
      std::cin >> R;
    }
  }
  catch (std::runtime_error& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  
  FrequencyMatrix F = R.toFrequencyMatrix(alpha, 2);
  std::cout << F;
}
