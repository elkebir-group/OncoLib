/*
 * tree2dotmain.cpp
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
  std::string filenameColorMap;
  std::string filenameVertexLabeling;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("c", "Color map", filenameColorMap)
    .refOption("l", "Vertex labeling", filenameVertexLabeling)
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
  
  StringToIntMap colorMap;
  if (filenameColorMap.empty())
  {
    colorMap = T.generateColorMap();
  }
  else
  {
    std::ifstream inColorMap(filenameColorMap.c_str());
    if (!inColorMap.good())
    {
      std::cerr << "Error: could not open '" << filenameColorMap << "' for reading" << std::endl;
      return 1;
    }
    
    if (!T.readColorMap(inColorMap, colorMap))
    {
      return 1;
    }
  }
  
  if (!filenameVertexLabeling.empty())
  {
    std::ifstream inVertexLabeling(filenameVertexLabeling.c_str());
    if (!inVertexLabeling.good())
    {
      std::cerr << "Error: could not open '" << filenameVertexLabeling << "' for reading" << std::endl;
      return 1;
    }
    
    StringNodeMap lPlus(T.tree());
    if (!T.readVertexLabeling(inVertexLabeling, T, lPlus))
    {
      return 1;
    }
    
    T.writeDOT(std::cout, lPlus, colorMap);
  }
  else
  {
    T.writeDOT(std::cout, colorMap);
  }
  
  return 0;
}

