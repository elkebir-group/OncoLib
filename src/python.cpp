/*
 *  python.cpp
 *
 *   Created on: 21-may-2018
 *       Author: M. El-Kebir
 */

#include <boost/scoped_array.hpp>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <iostream>
#include <fstream>
#include "utils.h"
#include "simulation.h"
#include "mutclonetree.h"
#include "readmatrix.h"
#include "generatemixture.h"

namespace p = boost::python;

template<typename T>
inline
std::vector< T > py_list_to_std_vector( const boost::python::object& iterable )
{
  return std::vector< T >( boost::python::stl_input_iterator< T >( iterable ),
                          boost::python::stl_input_iterator< T >( ) );
}

void simulate(std::string outputDirectory = ".",
              p::dict args = p::dict())
{
  std::string filenameColorMap;
  int seed = 0;
  double carryingCapacity = 5e4;
  double mutFreqThreshold = 0.05;
  double migrationRate = 1e-6;
  double mutationRate = 0.1;
  double driverRate = 2e-7;
  int maxNrAnatomicalSites = 3;
  int migrationPattern = 0;
  int nrTrials = -1;
  int sequencingDepth = 200;
  int nrSamplesPerAnatomicalSite = 2;
  int nrSamplesPrimary = 2;
  double sequencingErrorRate = 0;
  double purity = 1;
  bool verbose = false;
  
  py_list_to_std_vector<std::string>(args.keys());
  for (const std::string& key : py_list_to_std_vector<std::string>(args.keys()))
  {
    if (key == "seed")
    {
      seed = p::extract<int>(args.get(key));
    }
    else if (key == "carryingCapacity")
    {
      carryingCapacity = p::extract<double>(args.get(key));
    }
    else if (key == "mutFreqThreshold")
    {
      mutFreqThreshold = p::extract<double>(args.get(key));
    }
    else if (key == "migrationRate")
    {
      migrationRate = p::extract<double>(args.get(key));
    }
    else if (key == "mutationRate")
    {
      mutationRate = p::extract<double>(args.get(key));
    }
    else if (key == "driverRate")
    {
      driverRate = p::extract<double>(args.get(key));
    }
    else if (key == "maxNrAnatomicalSites")
    {
      maxNrAnatomicalSites = p::extract<int>(args.get(key));
    }
    else if (key == "migrationPattern")
    {
      migrationPattern = p::extract<int>(args.get(key));
    }
    else if (key == "nrTrials")
    {
      nrTrials = p::extract<int>(args.get(key));
    }
    else if (key == "sequencingDepth")
    {
      sequencingDepth = p::extract<int>(args.get(key));
    }
    else if (key == "nrSamplesPerAnatomicalSite")
    {
      nrSamplesPerAnatomicalSite = p::extract<int>(args.get(key));
    }
    else if (key == "nrSamplesPrimary")
    {
      nrSamplesPrimary = p::extract<int>(args.get(key));
    }
    else if (key == "sequencingErrorRate")
    {
      sequencingErrorRate = p::extract<double>(args.get(key));
    }
    else if (key == "purity")
    {
      purity = p::extract<double>(args.get(key));
    }
    else if (key == "filenameColorMap")
    {
      filenameColorMap = p::extract<std::string>(args.get(key));
    }
    else if (key == "verbose")
    {
      verbose = p::extract<bool>(args.get(key));
    }
  }
  
  Simulation::run(outputDirectory, filenameColorMap,
                  seed, carryingCapacity, mutFreqThreshold,
                  migrationRate, mutationRate, driverRate,
                  maxNrAnatomicalSites, migrationPattern, nrTrials,
                  sequencingDepth, nrSamplesPerAnatomicalSite, nrSamplesPrimary,
                  sequencingErrorRate, purity, verbose);
}

std::string visualizeMigrationGraph(std::string filenameInGraph,
                                    std::string filenameInColorMap = "")
{
  std::ifstream inG(filenameInGraph);
  if (!inG.good())
  {
    throw std::runtime_error("Error: could not open '"
                             + filenameInGraph + "' for reading");
  }
  
  MigrationGraph G;
  if (!G.read(inG))
  {
    throw std::runtime_error("Error: invalid format '"
                             + filenameInGraph + "'");
  }
  
  StringToIntMap colorMap;
  if (!filenameInColorMap.empty())
  {
    std::ifstream inColorMap(filenameInColorMap.c_str());
    if (!inColorMap.good())
    {
      throw std::runtime_error("Error: could not open '"
                               + filenameInColorMap + "' for reading");
    }
    if (!BaseTree::readColorMap(inColorMap, colorMap))
    {
      throw std::runtime_error("Error: error parsing '"
                               + filenameInColorMap + "'");
    }
  }
  else
  {
    colorMap = G.generateColorMap();
  }
  
  std::stringstream ss;
  
  G.writeDOT(ss, colorMap);
  
  return ss.str();
}

std::string visualizeCloneTree(std::string filenameInTree,
                               std::string filenameInColorMap = "",
                               std::string filenameInVertexLabeling = "")
{
  std::ifstream inT(filenameInTree.c_str());
  if (!inT.good())
  {
    throw std::runtime_error("Error: could not open '"
                             + filenameInTree + "' for reading");
  }
  
  MutCloneTree T;
  inT >> T;

  StringToIntMap colorMap;
  if (!filenameInColorMap.empty())
  {
    std::ifstream inColorMap(filenameInColorMap.c_str());
    if (!inColorMap.good())
    {
      throw std::runtime_error("Error: could not open '"
                               + filenameInColorMap + "' for reading");
    }
    if (!BaseTree::readColorMap(inColorMap, colorMap))
    {
      throw std::runtime_error("Error: error parsing '"
                               + filenameInColorMap + "'");
    }
  }
  else
  {
    colorMap = T.generateColorMap();
  }
  
  std::stringstream ss;
  
  if (filenameInVertexLabeling.empty())
  {
    T.writeDOT(ss, colorMap);
  }
  else
  {
    std::ifstream inVertexLabeling(filenameInVertexLabeling.c_str());
    if (!inVertexLabeling.good())
    {
      throw std::runtime_error("Error: could not open '"
                               + filenameInVertexLabeling
                               + "' for reading");
    }
    
    StringNodeMap lPlus(T.tree());
    if (!T.readVertexLabeling(inVertexLabeling, T, lPlus))
    {
      throw std::runtime_error("Error: error parsing '"
                               + filenameInVertexLabeling + "'");
    }
    
    T.writeDOT(ss, lPlus, colorMap);
  }
  
  return ss.str();
}

void mix(const std::string& filenameInTree,
         const std::string& filenameOutTree,
         int k,
         int seed = 0)
{
  std::ifstream inT(filenameInTree.c_str());
  if (!inT.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInTree + "' for reading");
  }
  
  MutCloneTree T;
  inT >> T;

  g_rng = std::mt19937(seed);
  MutCloneTree sampledT = GenerateMixture(T).generate(k);
  
  std::ofstream outT(filenameOutTree.c_str());
  outT << sampledT;
  outT.close();
}

void downSample(const std::string& filenameInR,
                const std::string& filenameOutR,
                int nrSamplesPerAnatomicalSite,
                int coverage,
                int seed = 0,
                double purity = 1,
                double sequencingErrorRate = 0,
                double fractionSNVs = 1)
{
  ReadMatrix R;
  std::ifstream inR(filenameInR.c_str());
  if (!inR.good())
  {
    throw std::runtime_error("Error: failed to open '" + filenameInR + "' for reading");
  }
  inR >> R;
  
  g_rng = std::mt19937(seed);
  ReadMatrix newR = R.downSample(nrSamplesPerAnatomicalSite,
                                 coverage,
                                 purity,
                                 sequencingErrorRate,
                                 fractionSNVs);
  
  std::ofstream outR(filenameOutR.c_str());
  outR << newR;
  outR.close();
}

void treeToFreqs(const std::string& filenameInTree,
                 const std::string& filenameOutFreqs)
{
  std::ifstream inT(filenameInTree.c_str());
  if (!inT.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInTree + "' for reading");
  }
  
  MutCloneTree T;
  inT >> T;
  
  std::ofstream outF(filenameOutFreqs.c_str());
  outF << T.getFrequencies();
  outF.close();
}

void readsToFreqs(const std::string& filenameInReads,
                  const std::string& filenameOutFreqs,
                  double alpha)
{
  std::ifstream inR(filenameInReads.c_str());
  if (!inR.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInReads + "' for reading");
  }
  
  ReadMatrix R;
  inR >> R;
  
  std::ofstream outF(filenameOutFreqs.c_str());
  outF << R.toFrequencyMatrix(alpha, 2);
  outF.close();
}

void sequence(const std::string& filenameInTree,
              const std::string& filenameOutReads,
              int sequencingDepth,
              int seed = 0,
              int ploidy = 2,
              double purity = 1.,
              double sequencingErrorRate = 0.)
{
  std::ifstream inT(filenameInTree.c_str());
  if (!inT.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInTree + "' for reading");
  }
  
  MutCloneTree T;
  inT >> T;
  
  g_rng = std::mt19937(seed);
  
  std::ofstream outR(filenameOutReads.c_str());
  outR << T.getReads(purity, sequencingDepth, sequencingErrorRate, ploidy);
  outR.close();
}

BOOST_PYTHON_FUNCTION_OVERLOADS(sequence_overloads, sequence, 3, 7);

BOOST_PYTHON_FUNCTION_OVERLOADS(downSample_overloads, downSample, 4, 8);

BOOST_PYTHON_FUNCTION_OVERLOADS(visualizeCloneTree_overloads, visualizeCloneTree, 1, 3);

BOOST_PYTHON_FUNCTION_OVERLOADS(visualizeMigrationGraph_overloads, visualizeMigrationGraph, 1, 2);

BOOST_PYTHON_FUNCTION_OVERLOADS(simulate_overloads, simulate, 0, 2);

BOOST_PYTHON_FUNCTION_OVERLOADS(mix_overloads, mix, 3, 4);

BOOST_PYTHON_MODULE(oncosim)
{
  p::def("sequence", sequence,
         sequence_overloads(p::args("filenameInTree",
                                    "filenameOutReads",
                                    "sequencingDepth",
                                    "seed=0",
                                    "ploidy=2",
                                    "purity=1.",
                                    "sequencingErrorRate=0."),
                            "Generate NGS reads from clone tree"));
  
  p::def("simulate", simulate, simulate_overloads(p::args("outputDirectory='.'",
                                                          "args={'filenameColorMap:'', "\
                                                          "'seed':0, "\
                                                          "'carryingCapacity':5e4, "\
                                                          "'mutFreqThreshold':0.05, "\
                                                          "'migrationRate':1e-6, "\
                                                          "'mutationRate':0.1, "\
                                                          "'driverRate':2e-7, "\
                                                          "'maxNrAnatomicalSites':3, "\
                                                          "'migrationPattern':0, "\
                                                          "'nrTrials':-1, "\
                                                          "'sequencingDepth':200, "\
                                                          "'nrSamplesPerAnatomicalSite':2, "\
                                                          "'nrSamplesPrimary':2, "\
                                                          "'sequencingErrorRate':0, "\
                                                          "'purity':1., "\
                                                          "'verbose':False}"),
                                                  "Simulate a metastatic tumor"));
  
  p::def("tree2dot", visualizeCloneTree,
         visualizeCloneTree_overloads(p::args("filenameInTree",
                                              "filenameInColorMap=''",
                                              "filenameInVertexLabeling=''"),
                                      "Visualize clone tree in DOT format"));
  
  p::def("graph2dot", visualizeMigrationGraph,
         visualizeMigrationGraph_overloads(p::args("filenameInTree",
                                                   "filenameInColorMap=''"),
                                      "Visualize clone tree in DOT format"));
  
  p::def("downsample", downSample,
         downSample_overloads(p::args("filenameInR",
                                      "filenameOutR",
                                      "nrSamplesPerAnatomicalSite",
                                      "coverage",
                                      "seed=0",
                                      "purity=1.",
                                      "sequencingErrorRate=0.",
                                      "fractionSNVs=1."),
                              "Down sample reads"));
  
  p::def("mix", mix, mix_overloads(p::args("filenameInTree",
                                           "filenameOutTree",
                                           "k",
                                           "seed=0"),
                                   "Generate mixture of clone tree leaves"));
  
  p::def("tree2freqs", treeToFreqs, "Extract mutation frequencies from clone tree");
  
  p::def("reads2freqs", readsToFreqs, "Extract mutation frequencies from reads");
}
