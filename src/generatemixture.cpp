/*
 *  generatemixture.cpp
 *
 *   Created on: 23-may-2018
 *       Author: M. El-Kebir
 */

#include "generatemixture.h"

GenerateMixture::GenerateMixture(const MutCloneTree& T)
  : _T(T)
{
}

MutCloneTree GenerateMixture::generate(int k,
                                       bool partition) const
{
  const int m = _T.getNrLocations();
  StringSet locationSet = _T.getLocations();
  StringVector locationVector(locationSet.begin(), locationSet.end());
  StringVector sampleVector, mutationVector;
  IntVector sampleToLocationVector;
  
  /// 1. Partition leaf set by locations and mutations
  StringToIntMap mutationMap;
  StringToNodeSetMap leavesByLocation;
  StringToNodeSetMap leavesByMutation;
  for (Node v : _T.leafSet())
  {
    const std::string& s = _T.l(v);
    
    leavesByLocation[s].insert(v);
    for (const std::string& j : _T.getMutations(v))
    {
      if (leavesByMutation.count(j) == 0)
      {
        leavesByMutation[j] = NodeSet();
      }
      leavesByMutation[j].insert(v);
      
      if (mutationMap.count(j) == 0)
      {
        int index = mutationMap.size();
        mutationMap[j] = index;
      }
    }
  }
  
  /// 2. Sample
  const int nrMutations = leavesByMutation.size();
  mutationVector = StringVector(nrMutations);
  for (const auto& pair : mutationMap)
  {
    const std::string& i = pair.first;
    int mapped_i = pair.second;

    mutationVector[mapped_i] = i;
  }
  
  DoubleVectorNodeMap sampleProportions(_T.tree(), DoubleVector(k, 0));
  
  if (partition)
  {
    IntNodeMap sampleIndexMap(_T.tree(), 0);
    std::uniform_int_distribution<> unif_int(0, k - 1);
    for (Node v : _T.leafSet())
    {
      sampleIndexMap[v] = unif_int(g_rng);
    }
    
    for (int s = 0; s < m; ++s)
    {
      const std::string& loc_s = locationVector[s];
      for (int p = 0; p < k; ++p)
      {
        double sum = 0;
        for (Node v : leavesByLocation[loc_s])
        {
          if (sampleIndexMap[v] == p)
          {
            sum += _T.getMixtureProportions(v)[0];
          }
        }
        
        for (Node v : leavesByLocation[loc_s])
        {
          if (sampleIndexMap[v] == p)
          {
            sampleProportions[v][p] = _T.getMixtureProportions(v)[0] / sum;
          }
        }
      }
    }
    return constructCloneTree(_T.leafSet(), sampleProportions);
  }
  else
  {
    NodeSet sampledLeaves;
    for (int s = 0; s < m; ++s)
    {
      const std::string& loc_s = locationVector[s];
      for (int p = 0; p < k; ++p)
      {
        // simulate draw from Dirichlet by drawing from gamma distributions
        // https://en.wikipedia.org/wiki/Dirichlet_distribution#Gamma_distribution
        
        DoubleNodeMap draw(_T.tree(), 0);
        double sum = 0;
        for (Node v : leavesByLocation[loc_s])
        {
          std::gamma_distribution<> gamma_dist(_T.getMixtureProportions(v)[0], 1);
          draw[v] = (gamma_dist(g_rng));
          sum += draw[v];
        }
        
        double new_sum = 0;
        for (Node v : leavesByLocation[loc_s])
        {
          if ((draw[v] / sum) >= 0.01)
          {
            sampleProportions[v][p] = draw[v];
            new_sum += draw[v];
            sampledLeaves.insert(v);
          }
          else
          {
            sampleProportions[v][p] = 0;
          }
        }
        
        for (Node v : leavesByLocation[loc_s])
        {
          sampleProportions[v][p] /= new_sum;
        }
      }
    }
    return constructCloneTree(sampledLeaves, sampleProportions);
  }
}

MutCloneTree GenerateMixture::constructCloneTree(const NodeSet& sampledLeaves,
                                                 const DoubleVectorNodeMap& sampleProportions) const
{
  Digraph newT;
  Node newRoot = lemon::INVALID;
  StringNodeMap locationStr(newT, "");
  StringNodeMap label(newT, "");
  
  BoolNodeMap filterNodes(_T.tree(), false);
  BoolArcMap filterArcs(_T.tree(), false);
  SubDigraph subT(_T.tree(), filterNodes, filterArcs);
  for (Node v : sampledLeaves)
  {
    filterNodes[v] = true;
    InArcIt a(_T.tree(), v);
    while (a != lemon::INVALID)
    {
      Node u = _T.tree().source(a);
      filterNodes[u] = true;
      filterArcs[a] = true;
      a = InArcIt(_T.tree(), u);
    }
  }
  
  /// The parameter should be a map, whose key type
  /// is the Node type of the destination digraph, while the value type is
  /// the Node type of the source digraph.
  NodeNodeMap ref(newT);
  
  StringSetNodeMap newMutationMap(newT);
  lemon::digraphCopy(subT, newT)
    .node(_T.root(), newRoot)
    .nodeMap(_T.getIdMap(), label)
    .nodeMap(_T.getMutationsNodeMap(), newMutationMap)
    .nodeMap(_T.getLeafLabeling(), locationStr)
    .nodeCrossRef(ref)
    .run();
  
  std::map<std::string, StringSet> labelToMutations;
  for (NodeIt v(newT); v != lemon::INVALID; ++v)
  {
    if (OutArcIt(newT, v) == lemon::INVALID)
    {
      labelToMutations[label[v]] = newMutationMap[v];
    }
  }
  
  // Contract degree one nodes
  bool changed = true;
  while (changed)
  {
    changed = false;
    for (NodeIt v(newT); v != lemon::INVALID; ++v)
    {
      if (v != newRoot && lemon::countOutArcs(newT, v) == 1)
      {
        Node u = newT.source(InArcIt(newT, v));
        Node w = newT.target(OutArcIt(newT, v));
        if (lemon::countOutArcs(newT, w) != 0)
        {
          if (labelToMutations.count(label[w]) == 1)
          {
            labelToMutations.erase(label[w]);
          }

          label[w] = label[v] + ";" + label[w];
          newMutationMap[w].insert(newMutationMap[v].begin(), newMutationMap[v].end());
          labelToMutations[label[w]] = newMutationMap[w];

          if (labelToMutations.count(label[v]) == 1)
          {
            labelToMutations.erase(label[v]);
          }

          newT.erase(v);

          newT.addArc(u, w);
          changed = true;
          break;
        }
      }
    }
  }
  
  MutCloneTree T(newT, newRoot, label, locationStr);
  for (NodeIt v(T.tree()); v != lemon::INVALID; ++v)
  {
    const std::string& label_v = T.label(v);
    if (T.isLeaf(v))
    {
      assert(labelToMutations.count(label_v) == 1);
      T.setMutations(v, labelToMutations[label_v]);
      Node vv = _T.getNodeByLabel(label_v);
      assert(vv != lemon::INVALID);
      T.setMixtureProportions(v, sampleProportions[vv]);
    }
  }
  
  return T;
}
