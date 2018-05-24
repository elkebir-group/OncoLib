/*
 * graph2dotmain.cpp
 *
 *  Created on: 24-may-2018
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include <lemon/arg_parser.h>
#include <fstream>
#include "migrationgraph.h"

int main(int argc, char** argv)
{
  std::string filenameColorMap;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("c", "Color map", filenameColorMap)
    .other("G", "Migration graph");
  ap.parse();
  
  MigrationGraph G;
  try
  {
    if (ap.files().size() >= 1)
    {
      std::string filenameGraph = ap.files()[0];
      std::ifstream inG(filenameGraph.c_str());
      if (!inG.good())
      {
        std::cerr << "Error: could not open '" << filenameGraph << "' for reading" << std::endl;
        return 1;
      }
      if (!G.read(inG))
      {
        return 1;
      }
    }
    else
    {
      if (!G.read(std::cin))
      {
        return 1;
      }
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
    colorMap = G.generateColorMap();
  }
  else
  {
    std::ifstream inColorMap(filenameColorMap.c_str());
    if (!inColorMap.good())
    {
      std::cerr << "Error: could not open '" << filenameColorMap << "' for reading" << std::endl;
      return 1;
    }
    
    if (!BaseTree::readColorMap(inColorMap, colorMap))
    {
      return 1;
    }
  }
  
  G.writeDOT(std::cout, colorMap);
  
  return 0;
}

