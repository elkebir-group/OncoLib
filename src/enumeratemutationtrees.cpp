/*
 * enumeratemutationtrees.cpp
 *
 *  Created on: 07-sep-2017
 *      Author: M. El-Kebir
 */

#include "enumeratemutationtrees.h"
#include "spruce/rootedcladisticnoisyenumeration.h"

EnumerateMutationTrees::EnumerateMutationTrees(const FrequencyMatrix& F)
  : _F(F)
  , _nrCharactersInTrees(0)
  , _canonicalStateTrees(_F.getNrCharacters(), StateTree({-1, 0}))
  , _F_lb(getFrequencyTensorLB())
  , _F_ub(getFrequencyTensorUB())
  , _G(_canonicalStateTrees, _F_lb, _F_ub)
{
  _G.init();
  _G.setLabels(_F_lb);
}

void EnumerateMutationTrees::enumerate(const std::string& outputDirectory,
                                       const int nrThreads,
                                       const int limit,
                                       const int timeLimit,
                                       const bool allCharacters,
                                       bool dryRun,
                                       std::ostream& out)
{
  int characterLowerBound = 1;
  if (allCharacters)
  {
    characterLowerBound = 0;
    for (int mult : _F.getCharacterMultiplicities())
    {
      characterLowerBound += mult;
    }
  }
  
  gm::RootedCladisticNoisyEnumeration enumerate(_G, limit, timeLimit,
                                                nrThreads,
                                                characterLowerBound,
                                                true,
                                                false,
                                                IntSet(),
                                                _F.getCharacterMultiplicities(),
                                                dryRun);
  
  enumerate.run();
  
  _nrCharactersInTrees = enumerate.objectiveValue();
  
  const gm::RootedCladisticNoisyEnumeration::ArcListList& solutions = enumerate.solutions();
  out << "# " << _nrCharactersInTrees << " out of "
      << _F.getMultiplicitySum()<< " mutations" << std::endl;
  out << (dryRun ? enumerate.getCounter() : solutions.size()) << " #trees" << std::endl;
  int index = 0;
  for (const auto& solution : solutions)
  {
    out << solution.size() - 1 << " #edges, tree " << index + 1 << std::endl;
    
    for (Arc a : solution)
    {
      Node u = _G.G().source(a);
      Node v = _G.G().target(a);
      
      if (u != _G.root())
      {
        int c = _G.nodeToCharStateList(v).front().first;
        int d = _G.nodeToCharStateList(u).front().first;
        out << _F.indexToCharacter(d) << " " << _F.indexToCharacter(c) << std::endl;
      }
    }
    ++index;
  }
  
  std::cerr << "Found " << (dryRun ? enumerate.getCounter() : solutions.size()) << " mutation trees with "
            << _nrCharactersInTrees << " out of "
            << _F.getMultiplicitySum()<< " mutations" << std::endl;
}


void EnumerateMutationTrees::enumerate(const std::string& outputDirectory,
                                       const int nrThreads,
                                       const int limit,
                                       const int timeLimit,
                                       const bool allCharacters,
                                       TreeVector& mutationTrees)
{
  int characterLowerBound = 1;
  if (allCharacters)
  {
    characterLowerBound = 0;
    for (int mult : _F.getCharacterMultiplicities())
    {
      characterLowerBound += mult;
    }
  }
  
  gm::RootedCladisticNoisyEnumeration enumerate(_G, limit, timeLimit,
                                                nrThreads,
                                                characterLowerBound,
                                                true,
                                                false,
                                                IntSet(),
                                                _F.getCharacterMultiplicities(),
                                                false);
  
  enumerate.run();
  
  _nrCharactersInTrees = enumerate.objectiveValue();
  
  const gm::RootedCladisticNoisyEnumeration::ArcListList& solutions = enumerate.solutions();
  for (const auto& solution : solutions)
  {
    Digraph newT;
    Node newRoot = lemon::INVALID;
    StringNodeMap idMap(newT, "");
    StringNodeMap l(newT, "");
    
    NodeVector nodeMap(solution.size(), lemon::INVALID);
    for (int i = 0; i < nodeMap.size(); ++i)
    {
      nodeMap[i] = newT.addNode();
    }
    
    for (Arc a : solution)
    {
      Node u = _G.G().source(a);
      Node v = _G.G().target(a);
      
      int c = _G.nodeToCharStateList(v).front().first;
      Node vv = nodeMap[c];
      idMap[vv] = _F.indexToCharacter(c);
      
      if (u == _G.root())
      {
        newRoot = vv;
        continue;
      }
      else
      {
        int d = _G.nodeToCharStateList(u).front().first;
        newT.addArc(nodeMap[d], vv);
      }
    }
    mutationTrees.push_back(CloneTree(newT, newRoot, idMap, l));
  }
  
  
  std::cerr << "Found " << solutions.size() << " mutation trees with "
            << _nrCharactersInTrees << " out of "
            << _F.getMultiplicitySum()<< " mutations" << std::endl;
}

void EnumerateMutationTrees::getFrequencyTensor(RealTensor& F_lb,
                                                RealTensor& F_ub) const
{
  const int n = _F.getNrCharacters();
  const int k = _F.getNrSamples();
  
  F_lb = RealTensor(2, k, n);
  F_ub = RealTensor(2, k, n);
  
  for (int p = 0; p < k; ++p)
  {
    F_lb.setRowLabel(p, _F.indexToSample(p));
    F_ub.setRowLabel(p, _F.indexToSample(p));
    for (int i = 0; i < n; ++i)
    {
      if (p == 0)
      {
        F_lb.setColLabel(i, _F.indexToCharacter(i));
        F_ub.setColLabel(i, _F.indexToCharacter(i));
      }
      F_lb.set(1, p, i, _F.min(p, i));
      F_ub.set(1, p, i, _F.max(p, i));
      F_lb.set(0, p, i, 1 - _F.max(p, i));
      F_ub.set(0, p, i, 1 - _F.min(p, i));
    }
  }
}

std::ostream& operator<<(std::ostream& out,
                         const EnumerateMutationTrees::TreeVector& trees)
{
  out << trees.size() << " #trees" << std::endl;
  int idx = 0;
  for (const CloneTree& T : trees)
  {
    out << lemon::countArcs(T.tree()) << " #edges, tree " << idx + 1 << std::endl;
    T.write(out);
    ++idx;
  }
  return out;
}

std::istream& operator>>(std::istream& in,
                         EnumerateMutationTrees::TreeVector& trees)
{
  int nrTrees = -1;
  
  std::string line;
  getline(in, line);
  
  while (line.empty() || line[0] == '#')
  {
    getline(in, line);
  }
  
  std::stringstream ss(line);
  ss >> nrTrees;
  
  if (nrTrees < 0)
  {
    throw std::runtime_error("Error: number of trees should be nonnegative");
  }
  
  for (int i = 0; i < nrTrees; ++i)
  {
    int nrEdges = -1;
    
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
    
    trees.push_back(T);
  }
  
  return in;
}
