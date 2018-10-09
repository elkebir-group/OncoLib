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
#include <boost/python/stl_iterator.hpp>
#include <armadillo>
#include <iostream>
#include <fstream>
#include "utils.h"
#include "simulation.h"
#include "mutclonetree.h"
#include "readmatrix.h"
#include "generatemixture.h"
#include "enumeratemutationtrees.h"

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
  int desiredNrMutationClusters = -1;
  int sequencingDepth = 200;
  int nrSamplesPerAnatomicalSite = 2;
  int nrSamplesPrimary = 2;
  double sequencingErrorRate = 0;
  double purity = 1;
  bool verbose = false;
  
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
    else if (key == "desiredNrMutationClusters")
    {
      desiredNrMutationClusters = p::extract<int>(args.get(key));
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
                  maxNrAnatomicalSites, migrationPattern,
                  nrTrials, desiredNrMutationClusters,
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
         bool partition,
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
  MutCloneTree sampledT = GenerateMixture(T).generate(k, partition);
  
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
  inR.close();
  
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
  outR << T.getReads(purity, sequencingDepth, sequencingErrorRate, ploidy, true, true);
  outR.close();
}

void precluster(const std::string& filenameInFreqs,
                const std::string& filenameInClustering,
                const std::string& filenameOutFreqs)
{
  std::ifstream inF(filenameInFreqs.c_str());
  if (!inF.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInFreqs + "' for reading");
  }
  
  FrequencyMatrix F;
  inF >> F;
  inF.close();
  
  std::ifstream inC(filenameInClustering.c_str());
  if (!inC.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInClustering + "' for reading");
  }
  
  IntMatrix clustering = F.parseClustering(inC);
  inC.close();
  
  FrequencyMatrix newF = F.cluster(clustering);
  
  std::ofstream outF(filenameOutFreqs);
  outF << newF;
  outF.close();
}

void enumerate(const std::string& filenameInFreqs,
               const std::string& filenameOutTrees,
               bool spanning,
               int nrThreads,
               int timeLimit,
               bool dryRun)
{
  std::ifstream inF(filenameInFreqs.c_str());
  if (!inF.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInFreqs + "' for reading");
  }
  
  FrequencyMatrix F;
  inF >> F;
  inF.close();
  
  VerbosityLevel oldVerbosity = g_verbosity;
  g_verbosity = VERBOSE_NONE;

  std::ofstream outT(filenameOutTrees.c_str());
  
  EnumerateMutationTrees enumerate(F);
  enumerate.enumerate("",
                      nrThreads,
                      -1,
                      timeLimit,
                      spanning,
                      dryRun,
                      outT);
  
  g_verbosity = oldVerbosity;
  
  outT.close();
}

double countSpanningTrees(const std::string& filenameInFreqs,
                          int root = 0)
{
  std::ifstream inF(filenameInFreqs.c_str());
  if (!inF.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInFreqs + "' for reading");
  }
  
  FrequencyMatrix F;
  inF >> F;
  inF.close();
  
  EnumerateMutationTrees mutT(F);
  
  const Digraph& G = mutT.getAncestryGraph().G();
  lemon::DynArcLookUp<Digraph> arcLookUp(G);
  
  arma::mat L(F.getNrCharacters() - 1, F.getNrCharacters() - 1);
  
  for (int i = 0; i < F.getNrCharacters(); ++i)
  {
    for (int j = 0; j < F.getNrCharacters(); ++j)
    {
      if (i == root || j == root)
      {
        continue;
      }
      
      int ii = i > root ? i - 1 : i;
      int jj = j > root ? j - 1 : j;
      
      if (i == j)
      {
        Node v_j = mutT.getAncestryGraph().charStateToNode(j, 1);
        L(ii, jj) = lemon::countInArcs(G, v_j) - 1;
      }
      else
      {
        Node v_i = mutT.getAncestryGraph().charStateToNode(i, 1);
        Node v_j = mutT.getAncestryGraph().charStateToNode(j, 1);
        
        if (arcLookUp(v_i, v_j) != lemon::INVALID)
        {
          L(ii, jj) = -1.;
        }
        else
        {
          L(ii, jj) = 0.;
        }
      }
    }
  }
  
//  std::cout << L << std::endl;
  
  return arma::det(L);
}

CloneTree parseNextTree(std::istream& in)
{
  int nrEdges = -1;
  std::string line;
  std::stringstream ss;
  
  getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> nrEdges;
  
  if (nrEdges < 0)
  {
    throw std::runtime_error("Error: number of edges should be nonnegative");
  }
  
  ss.clear();
  ss.str("");
  for (int j = 0; j < nrEdges; ++j)
  {
    getline(in, line);
    ss << line << std::endl;
  }
  
  CloneTree T;
  if (!T.read(ss))
    throw std::runtime_error("");
  
  return T;
}

void filterSCS(const std::string& filenameInInferredTrees,
               const std::string& filenameInTrueTree,
               const std::string& filenameOutTrees,
               int nrClones,
               int seed)
{
  // 0. parse true clone tree
  std::ifstream inTrueTree(filenameInTrueTree.c_str());
  if (!inTrueTree.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInTrueTree + "' for reading");
  }
  
  MutCloneTree trueT;
  inTrueTree >> trueT;
  inTrueTree.close();
  
  // 1. collect clones
  std::vector<StringVector> clones;
  for (Node leaf : trueT.leafSet())
  {
    clones.push_back(StringVector());
    Node u = leaf;
    while ((u = trueT.parent(u)) != trueT.root())
    {
      clones.back().push_back(trueT.label(u));
    }
  }
  
  g_rng = std::mt19937(seed);
  std::shuffle(clones.begin(), clones.end(), g_rng);
  while (clones.size() > nrClones)
  {
    clones.pop_back();
  }
  
  // 2. parse enumerated trees
  std::ifstream inTrees(filenameInInferredTrees.c_str());
  if (!inTrees.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInInferredTrees + "' for reading");
  }
  
  int nrTrees = -1;
  std::string line;
  while (line.empty() || line[0] == '#')
  {
    getline(inTrees, line);
  }
  std::stringstream ss(line);
  ss >> nrTrees;
  
  if (nrTrees < 0)
  {
    throw std::runtime_error("Error: number of trees should be nonnegative");
  }
  
  std::ofstream outTrees(filenameOutTrees);
  for (const StringVector& clone : clones)
  {
    outTrees << "#";
    for (const std::string& mut : clone)
    {
      outTrees << " " << mut;
    }
    outTrees << std::endl;
  }
  
  for (int i = 0; i < nrTrees; ++i)
  {
    StringPairSet diff;
    CloneTree T = parseNextTree(inTrees);
    
    bool ok = true;
    for (const StringVector& clone : clones)
    {
      Node u = T.getNodeByLabel(clone[0]);
      int idx = 0;
      while (u != lemon::INVALID && ok)
      {
        ok &= (idx < clone.size()) && (T.label(u) == clone[idx]);
        u = T.parent(u);
        ++idx;
      }
    }
    
    if (ok)
    {
      outTrees << lemon::countArcs(T.tree()) << " #edges" << std::endl;
      T.write(outTrees);
    }
  }
  outTrees.close();
}

void filterLR(const std::string& filenameInInferredTrees,
              const std::string& filenameInTrueTree,
              const std::string& filenameOutTrees,
              int nrPairs,
              int seed)
{
  // 0. parse true clone tree
  std::ifstream inTrueTree(filenameInTrueTree.c_str());
  if (!inTrueTree.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInTrueTree + "' for reading");
  }
  
  MutCloneTree trueT;
  inTrueTree >> trueT;
  inTrueTree.close();
  
  // 1. collect clones
  StringPairSet ancestralPairSet = trueT.getAncestralPairsNoRoot();
  std::vector<StringPair> vec(ancestralPairSet.begin(), ancestralPairSet.end());
  
  g_rng = std::mt19937(seed);
  std::shuffle(vec.begin(), vec.end(), g_rng);
  while (vec.size() > nrPairs)
  {
    vec.pop_back();
  }
  ancestralPairSet = StringPairSet(vec.begin(), vec.end());
  
  // 2. parse enumerated trees
  std::ifstream inTrees(filenameInInferredTrees.c_str());
  if (!inTrees.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInInferredTrees + "' for reading");
  }
  
  int nrTrees = -1;
  std::string line;
  while (line.empty() || line[0] == '#')
  {
    getline(inTrees, line);
  }
  std::stringstream ss(line);
  ss >> nrTrees;
  
  if (nrTrees < 0)
  {
    throw std::runtime_error("Error: number of trees should be nonnegative");
  }
  
  std::ofstream outTrees(filenameOutTrees);
  for (const StringPair& pair : ancestralPairSet)
  {
    outTrees << "# " << pair.first << " " << pair.second << std::endl;
  }
  
  for (int i = 0; i < nrTrees; ++i)
  {
    StringPairSet diff;
    CloneTree T = parseNextTree(inTrees);
    StringPairSet ancestralPairSetT = T.getAncestralPairs();
    
    std::set_difference(ancestralPairSet.begin(), ancestralPairSet.end(),
                        ancestralPairSetT.begin(), ancestralPairSetT.end(),
                        std::inserter(diff, diff.begin()));
    
    if (diff.empty())
    {
      outTrees << lemon::countArcs(T.tree()) << " #edges" << std::endl;
      T.write(outTrees);
    }
  }
  outTrees.close();
}

void computeRecall(const std::string& filenameInInferredTrees,
                   const std::string& filenameInTrueTree,
                   const std::string& filenameOutRecall)
{
  std::ifstream inTrueTree(filenameInTrueTree.c_str());
  if (!inTrueTree.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInTrueTree + "' for reading");
  }
  
  CloneTree trueT;
  if (!trueT.read(inTrueTree))
  {
    throw std::runtime_error("Parse error: '" + filenameInTrueTree + "'");
  }
  inTrueTree.close();
  
  StringPairSet trueParental = trueT.getParentalPairs();
  StringPairSet trueAncestral = trueT.getAncestralPairs();
  StringPairSet trueIncomparable = trueT.getIncomparablePairs();
  
  std::ifstream inTrees(filenameInInferredTrees.c_str());
  if (!inTrees.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInInferredTrees + "' for reading");
  }
  
  int nrTrees = -1;
  std::string line;
  while (line.empty() || line[0] == '#')
  {
    getline(inTrees, line);
  }
  
  std::stringstream ss(line);
  ss >> nrTrees;
  if (nrTrees < 0)
  {
    throw std::runtime_error("Error: number of trees should be nonnegative");
  }
  
  std::ofstream outRecall(filenameOutRecall);
  outRecall << "solution\tancestral\tparental\tincomparable" << std::endl;
  for (int treeIdx = 0; treeIdx < nrTrees; ++treeIdx)
  {
    CloneTree T = parseNextTree(inTrees);
    StringPairSet inferredParental = T.getParentalPairs();
    StringPairSet inferredAncestral = T.getAncestralPairs();
    StringPairSet inferredIncomparable = T.getIncomparablePairs();
    
    outRecall << treeIdx << "\t"
              << BaseTree::recall(inferredAncestral, trueAncestral)
              << "\t"
              << BaseTree::recall(inferredParental, trueParental)
              << "\t"
              << BaseTree::recall(inferredIncomparable, trueIncomparable)
              << std::endl;
  }
  inTrees.close();
  outRecall.close();
}

double getFractionOfIncomparablePairs(const std::string& filenameInFreqs)
{
  std::ifstream inF(filenameInFreqs.c_str());
  if (!inF.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInFreqs + "' for reading");
  }
  
  FrequencyMatrix F;
  inF >> F;
  inF.close();
  
  return EnumerateMutationTrees(F).getAncestryGraph().fracOfIncomparablePairs();
}

void summarize(const std::string& filenameInAllTrees,
               const std::string& filenameInSampledTree,
               const std::string& filenameOutSummary)
{
  // 1. Parse sampled trees
  std::map<StringPairSet, int> sampledTreeHistogram;
  
  std::ifstream inTrees(filenameInSampledTree.c_str());
  if (!inTrees.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInSampledTree + "' for reading");
  }
  
  int nrTrees = -1;
  std::string line;
  while (line.empty() || line[0] == '#')
  {
    getline(inTrees, line);
  }
  
  std::stringstream ss(line);
  ss >> nrTrees;
  if (nrTrees < 0)
  {
    throw std::runtime_error("Error: number of trees should be nonnegative");
  }
  
  for (int treeIdx = 0; treeIdx < nrTrees; ++treeIdx)
  {
    CloneTree T = parseNextTree(inTrees);
    StringPairSet edges;
    for (ArcIt a(T.tree()); a != lemon::INVALID; ++a)
    {
      Node u = T.tree().source(a);
      Node v = T.tree().target(a);
      
      edges.insert(StringPair(T.label(u), T.label(v)));
    }
    if (sampledTreeHistogram.count(edges) == 0)
    {
      sampledTreeHistogram[edges] = 1;
    }
    else
    {
      ++sampledTreeHistogram[edges];
    }
  }
  inTrees.close();
  
  // 2. Parse all trees
  std::map<StringPairSet, int> allTreeIndex;
  
  inTrees = std::ifstream(filenameInAllTrees.c_str());
  if (!inTrees.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInAllTrees + "' for reading");
  }
  
  nrTrees = -1;
  line.clear();
  while (line.empty() || line[0] == '#')
  {
    getline(inTrees, line);
  }
  
  ss.clear();
  ss.str(line);
  ss >> nrTrees;
  if (nrTrees < 0)
  {
    throw std::runtime_error("Error: number of trees should be nonnegative");
  }
  
  for (int treeIdx = 0; treeIdx < nrTrees; ++treeIdx)
  {
    CloneTree T = parseNextTree(inTrees);
    StringPairSet edges;
    for (ArcIt a(T.tree()); a != lemon::INVALID; ++a)
    {
      Node u = T.tree().source(a);
      Node v = T.tree().target(a);
      
      edges.insert(StringPair(T.label(u), T.label(v)));
    }
    if (allTreeIndex.count(edges) > 0)
    {
      throw std::runtime_error("Error: duplicate tree");
    }
    
    allTreeIndex[edges] = treeIdx;
  }
  inTrees.close();
  
  std::ofstream outSum(filenameOutSummary);
  int incorrect = 0;
  outSum << "index\tcount" << std::endl;
  for (const auto& kv : sampledTreeHistogram)
  {
    if (allTreeIndex.count(kv.first) == 0)
    {
      outSum << --incorrect << "\t" << kv.second << std::endl;
    }
    else
    {
      outSum << allTreeIndex[kv.first] << "\t" << kv.second << std::endl;
    }
  }
  for (const auto& kv : allTreeIndex)
  {
    if (sampledTreeHistogram.count(kv.first) == 0)
    {
      outSum << kv.second << "\t" << 0 << std::endl;
    }
  }
  outSum.close();
}

int rejectionSample(const std::string& filenameInFreqs,
                    const std::string& filenameOutTrees,
                    int count,
                    int seed = 0)
{
  g_rng = std::mt19937(seed);
  
  std::ifstream inF(filenameInFreqs.c_str());
  if (!inF.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInFreqs + "' for reading");
  }
  
  FrequencyMatrix F;
  inF >> F;
  inF.close();
  
  const int n = F.getNrCharacters();
  const int k = F.getNrSamples();
  
  RealTensor FF(2, k, n);
  for (int p = 0; p < k; ++p)
  {
    FF.setRowLabel(p, F.indexToSample(p));
    for (int i = 0; i < n; ++i)
    {
      if (p == 0)
      {
        FF.setColLabel(i, F.indexToCharacter(i));
      }
      FF.set(1, p, i, F.min(p, i));
      FF.set(0, p, i, 1 - F.max(p, i));
    }
  }
  
  StateTreeVector S(F.getNrCharacters(), StateTree({-1, 0}));
  
  gm::RootedCladisticAncestryGraph G(FF, S);
  G.init();
  
  std::ofstream outTrees(filenameOutTrees);
  outTrees << count << " #trees" << std::endl;
  StringPairSet tree;
  int sampleCount = 0;
  while (count > 0)
  {
    ++sampleCount;
    if (G.sample(tree))
    {
      outTrees << tree.size() << " #edges" << std::endl;
      for (const StringPair& edge : tree)
      {
        outTrees << edge.first << " " << edge.second << std::endl;
      }
      --count;
    }
  }
  
  return sampleCount;
}

std::string visualizeEnumeratedCloneTree(const std::string& filenameInFreqs,
                                         const std::string& filenameInTrees,
                                         int index,
                                         bool showCharacterLabels = false,
                                         int offset = 1,
                                         bool showFrequencies = false)
{
  std::ifstream inF(filenameInFreqs.c_str());
  if (!inF.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInFreqs + "' for reading");
  }
  
  FrequencyMatrix F;
  inF >> F;
  inF.close();
  
  std::ifstream inTrees(filenameInTrees.c_str());
  if (!inTrees.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInTrees + "' for reading");
  }
  
  int nrTrees = -1;
  std::string line;
  while (line.empty() || line[0] == '#')
  {
    getline(inTrees, line);
  }
  
  std::stringstream ss(line);
  ss >> nrTrees;
  if (nrTrees < 0)
  {
    throw std::runtime_error("Error: number of trees should be nonnegative");
  }
  
  if (!(0 <= index && index < nrTrees))
  {
    throw std::runtime_error("Error: invalid index");
  }

  std::stringstream ssOut;
  for (int treeIdx = 0; treeIdx < nrTrees; ++treeIdx)
  {
    CloneTree T = parseNextTree(inTrees);
    if (index == treeIdx)
    {
      ssOut << "digraph T_" << index << " {" << std::endl;
      
      for (NodeIt v(T.tree()); v != lemon::INVALID; ++v)
      {
        ssOut << "  " << T.tree().id(v) << " [label=\"";
        
        const std::string& label_v = T.label(v);
        if (showCharacterLabels)
          ssOut << label_v;
        else
          ssOut << F.characterToIndex(label_v) + offset;
        
        if (showFrequencies)
        {
          for (int p = 0; p < F.getNrSamples(); ++p)
          {
            ssOut << "\\n" << F.min(p, F.characterToIndex(label_v)) << " " << F.max(p, F.characterToIndex(label_v));
          }
        }
        
        ssOut << "\"]" << std::endl;
      }
      
      for (ArcIt a(T.tree()); a != lemon::INVALID; ++a)
      {
        Node u = T.tree().source(a);
        Node v = T.tree().target(a);
        
        ssOut << "  " << T.tree().id(u) << " -> " << T.tree().id(v) << std::endl;
      }
      
      ssOut << "}" << std::endl;
      break;
    }
  }
  inTrees.close();
  
  return ssOut.str();
}

int identifySolution(const std::string& filenameInTrees,
                     const std::string& filenameInTree)
{
  std::ifstream inTrueTree(filenameInTree.c_str());
  if (!inTrueTree.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInTree + "' for reading");
  }
  
  CloneTree trueT;
  if (!trueT.read(inTrueTree))
  {
    throw std::runtime_error("Parse error: '" + filenameInTree + "'");
  }
  inTrueTree.close();
  StringPairSet edgesTrueT = trueT.getEdgeSet();
  
  std::ifstream inTrees(filenameInTrees.c_str());
  if (!inTrees.good())
  {
    throw std::runtime_error("Error: could not open '" + filenameInTrees + "' for reading");
  }
  
  int nrTrees = -1;
  std::string line;
  while (line.empty() || line[0] == '#')
  {
    getline(inTrees, line);
  }
  
  std::stringstream ss(line);
  ss >> nrTrees;
  if (nrTrees < 0)
  {
    throw std::runtime_error("Error: number of trees should be nonnegative");
  }
  
  int res = -1;
  for (int treeIdx = 0; treeIdx < nrTrees; ++treeIdx)
  {
    CloneTree T = parseNextTree(inTrees);
    if (T.getEdgeSet() == edgesTrueT)
    {
      res = treeIdx;
      break;
    }
  }
  inTrees.close();
  
  return res;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(sequence_overloads, sequence, 3, 7);

BOOST_PYTHON_FUNCTION_OVERLOADS(downSample_overloads, downSample, 4, 8);

BOOST_PYTHON_FUNCTION_OVERLOADS(visualizeCloneTree_overloads, visualizeCloneTree, 1, 3);

BOOST_PYTHON_FUNCTION_OVERLOADS(visualizeEnumeratedCloneTree_overloads, visualizeEnumeratedCloneTree, 3, 6);

BOOST_PYTHON_FUNCTION_OVERLOADS(visualizeMigrationGraph_overloads, visualizeMigrationGraph, 1, 2);

BOOST_PYTHON_FUNCTION_OVERLOADS(simulate_overloads, simulate, 0, 2);

BOOST_PYTHON_FUNCTION_OVERLOADS(countSpanningTrees_overloads, countSpanningTrees, 1, 2);

BOOST_PYTHON_FUNCTION_OVERLOADS(mix_overloads, mix, 4, 5);

BOOST_PYTHON_FUNCTION_OVERLOADS(rejectionSample_overloads, rejectionSample, 3, 4);

BOOST_PYTHON_MODULE(oncolib)
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
                                                          "'desiredNrMutationClusters':-1, "
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
  
  p::def("enumeratedtree2dot", visualizeEnumeratedCloneTree,
         visualizeEnumeratedCloneTree_overloads(p::args("filenameInFreqs",
                                                        "filenameInTrees",
                                                        "index",
                                                        "showCharacterLabels=False",
                                                        "offset=1",
                                                        "showFrequencies=False"),
                                                "Visualize enumerated clone tree in DOT format"));
  
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
                                           "partition",
                                           "seed=0"),
                                   "Generate mixture of clone tree leaves"));
  
  p::def("tree2freqs", treeToFreqs, "Extract mutation frequencies from clone tree");
  
  p::def("reads2freqs", readsToFreqs, "Extract mutation frequencies from reads");
  
  p::def("precluster", precluster, "Cluster mutation frequencies given clustering");
  
  p::def("enumerate", enumerate, "Enumerate mutation trees");
  
  p::def("countSpanningTrees", countSpanningTrees,
         countSpanningTrees_overloads(p::args("filenameInF",
                                              "root=0"),
                                      "Compute number of spanning arborescence in ancestry graph"));
  
  p::def("computeRecall", computeRecall, "Compute recall");
  
  p::def("getFractionOfIncomparablePairs", getFractionOfIncomparablePairs,
         "Compute fraction of incomparable pairs");
  
  p::def("filterSCS", filterSCS, "Include SCS information");
  
  p::def("filterLR", filterLR, "Include long read information");
  
  p::def("summarize", summarize, "Compute histogram of sampled trees w.r.t. all trees");
  
  p::def("identifySolution", identifySolution, "Determine whether specified tree occurs in specified trees");
  
  p::def("rejectionSample", rejectionSample,
         rejectionSample_overloads(p::args("filenameInFreqs",
                                           "filenameOutTrees",
                                           "count",
                                           "seed=0"),
                                   "Rejection sampling of spanning trees satisfying SC"));
}
