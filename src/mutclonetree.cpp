/*
 * mutclonetree.cpp
 *
 *  Created on: 21-may-2018
 *      Author: M. El-Kebir
 */

#include "mutclonetree.h"
#include <iomanip>

MutCloneTree::MutCloneTree()
  : CloneTree()
  , _U(_tree)
  , _b(_tree)
{
}

MutCloneTree::MutCloneTree(const MutCloneTree& other)
  : CloneTree(other)
  , _U(_tree)
  , _b(_tree)
{
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    const std::string& lbl = label(u);
    Node other_u = other.getNodeByLabel(lbl);
    
    _U[u] = other._U[other_u];
    _b[u] = other._b[other_u];
  }
}

MutCloneTree& MutCloneTree::operator =(const MutCloneTree& other)
{
  if (this != &other)
  {
    lemon::digraphCopy(other._tree, _tree)
      .node(other._root, _root)
      .nodeMap(other._nodeToId, _nodeToId)
      .run();
    
    for (NodeIt u(_tree); u != lemon::INVALID; ++u)
    {
      const std::string& str = _nodeToId[u];
      _idToNode[str] = u;
    }
    
    init();
    
    for (NodeIt u(_tree); u != lemon::INVALID; ++u)
    {
      const std::string& lbl = label(u);
      Node other_u = other.getNodeByLabel(lbl);
      
      _l[u] = other._l[other_u];
      _U[u] = other._U[other_u];
      _b[u] = other._b[other_u];
    }
  }
  
  return *this;
}

MutCloneTree::MutCloneTree(const Digraph& T,
                           Node root,
                           const StringNodeMap& label,
                           const StringNodeMap& l)
  : CloneTree(T, root, label, l)
  , _U(_tree)
  , _b(_tree)
{
  for (NodeIt v(T); v != lemon::INVALID; ++v)
  {
    if (OutArcIt(T, v) == lemon::INVALID)
    {
      const std::string& label_v = label[v];
      _l[getNodeByLabel(label_v)] = l[v];
    }
  }
}

bool MutCloneTree::readLeafLabeling(std::istream& in)
{
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    _l[u] = "";
  }
  
  while (in.good())
  {
    std::string line;
    getline(in, line);
    
    if (line.empty())
      continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t "));
    
    if (s.size() < 3)
    {
      std::cerr << "Error: invalid format" << std::endl;
      return false;
    }
    
    std::string label_u = s[0];
    std::string label_s = s[1];
    
    if (_idToNode.count(label_u) == 0)
    {
      std::cerr << "Error: clone-tree vertex with label '" << label_u << "' does not exist" << std::endl;
      return false;
    }
    
    Node u = _idToNode[label_u];
    _l[u] = label_s;
    int nrSamples = boost::lexical_cast<int>(s[2]);
    if (s.size() < 3 + nrSamples)
    {
      std::cerr << "Error: invalid number of samples" << std::endl;
      return false;
    }
    
    _U[u] = DoubleVector(nrSamples, 0);
    for (int p = 0; p < nrSamples; ++p)
    {
      _U[u][p] = boost::lexical_cast<double>(s[3 + p]);
    }

    int nrCharacters = s.size() - 3 - nrSamples;
    _b[u].clear();
    for (int c = 0; c < nrCharacters; ++c)
    {
      _b[u].insert(boost::lexical_cast<int>(s[3 + nrSamples + c]));
    }
  }
  
  for (Node u : _leafSet)
  {
    if (_l[u].empty())
    {
      std::cerr << "Error: leaf '" << label(u) << "' left unlabeled" << std::endl;
      return false;
    }
  }
  
  return true;
}

void MutCloneTree::writeLeafLabeling(std::ostream& out) const
{
  for (Node u : _leafSet)
  {
    out << _nodeToId[u] << " " << _l[u] << " " << _U[u].size();
    for (double proportion : _U[u])
    {
      out << " " << proportion;
    }
    for (int c : _b[u])
    {
      out << " " << c;
    }
    out << std::endl;
  }
}

void MutCloneTree::writeDOT(std::ostream& out,
                            const StringToIntMap& colorMap) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(l(u))->second << ",label=\""
          << _nodeToId[u] << "\\n"
          << _l[u] << "\\n";
      
      bool first = true;
      for (double usage : _U[u])
      {
        if (first)
          first = false;
        else
          out << " ";
        out << std::setprecision(2) << 100 * usage << "%";
      }
      out << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [label=\"" << _nodeToId[u] << "\"]" << std::endl;
    }
  }
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    out << "\t" << _tree.id(_tree.source(a)) << " -> " << _tree.id(_tree.target(a)) << std::endl;
  }
  
  out << "}" << std::endl;
}

void MutCloneTree::writeDOT(std::ostream& out,
                            const StringNodeMap& lPlus,
                            const StringToIntMap& colorMap) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(l(u))->second << ",label=\""
          << _nodeToId[u] << "\\n"
          << _l[u] << "\\n";
      
      bool first = true;
      for (double usage : _U[u])
      {
        if (first)
          first = false;
        else
          out << " ";
        out << std::setprecision(2) << 100 * usage << "%";
      }
      out << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second << ",label=\"" << _nodeToId[u] << "\"]" << std::endl;
    }
  }
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    Node u = _tree.source(a);
    Node v = _tree.target(a);
    
    const std::string& s_u = lPlus[u];
    const std::string& s_v = lPlus[v];
    
    out << "\t" << _tree.id(u) << " -> " << _tree.id(v);
    if (s_u == s_v)
    {
      out << " [penwidth=3,colorscheme=set19,color=" << colorMap.find(s_u)->second << "]";
    }
    else
    {
      out << " [penwidth=3,colorscheme=set19,color=\"" << colorMap.find(s_u)->second << ";0.5:" << colorMap.find(s_v)->second << "\"]";
    }
    
    out << std::endl;
  }
  
  out << "}" << std::endl;
}

FrequencyMatrix MutCloneTree::getFrequencies() const
{
  char buf[1024];
  
  StringVector indexToLocation;
  StringVector indexToSample;
  StringVector indexToCharacter;
  
  IntSet characters;
  for (Node u : leafSet())
  {
    characters.insert(getMutations(u).begin(), getMutations(u).end());
  }
  
  for (int i : characters)
  {
    snprintf(buf, 1024, "%d", i);
    indexToCharacter.push_back(buf);
  }
  
  StringSet locationSet = getLocations();
  for (const std::string& location : locationSet)
  {
    indexToLocation.push_back(location);
  }
  
  IntVector sampleIndexToLocationIndex;
  IntVector sampleIndexToOrgSampleIndex;
  for (int s = 0; s < locationSet.size(); ++s)
  {
    const std::string& location = indexToLocation[s];
    for (Node u : leafSet())
    {
      if (l(u) == location)
      {
        const int nrSamples_s = _U[u].size();
        for (int p = 0; p < nrSamples_s; ++p)
        {
          snprintf(buf, 1024, "%s_%d", location.c_str(), p);
          indexToSample.push_back(buf);
          sampleIndexToLocationIndex.push_back(s);
          sampleIndexToOrgSampleIndex.push_back(p);
        }
        break;
      }
    }
  }
  
  FrequencyMatrix F(indexToLocation, indexToSample, indexToCharacter,
                    sampleIndexToLocationIndex);
  
  for (int i : characters)
  {
    snprintf(buf, 1024, "%d", i);
    const int c = F.characterToIndex(buf);
    
    for (Node u : leafSet())
    {
      if (hasMutation(u, i))
      {
        const std::string& location = l(u);
        const int s = F.anatomicalSiteToIndex(location);
        
        for (int p : F.anatomicalSiteIndexToSampleIndices(s))
        {
          double f = F.min(p, c) + _U[u][sampleIndexToOrgSampleIndex[p]];
          f = std::min(1., f);
          F.set(p, c, f, f);
        }
      }
    }
  }
  
  return F;
}

ReadMatrix MutCloneTree::getReads(double purity,
                                  int sequencingDepth,
                                  double sequencingErrorRate,
                                  int ploidy) const
{
  FrequencyMatrix F = getFrequencies();
  ReadMatrix R(F.getIndexToAnatomicalSites(),
               F.getIndexToSamples(),
               F.getIndexToCharacters(),
               F.getSampleIndexToAnatomicalSiteIndex());
  
  std::poisson_distribution<> poisson(sequencingDepth);
  for (int s = 0; s < F.getNrAnatomicalSites(); ++s)
  {
    for (int p : F.anatomicalSiteIndexToSampleIndices(s))
    {
      for (int c = 0; c < F.getNrCharacters(); ++c)
      {
        double f = F.min(p, c) * purity / ploidy;
        
        int coverage = poisson(g_rng);
        std::binomial_distribution<> binom(coverage, f);
        
        int org_var = binom(g_rng);
        int org_ref = coverage - org_var;
        
        if (g_tol.nonZero(sequencingErrorRate))
        {
          std::binomial_distribution<> binom_noise_var(org_var,
                                                       sequencingErrorRate);
          std::binomial_distribution<> binom_noise_ref(org_ref,
                                                       sequencingErrorRate);
          
          int flips_var =  binom_noise_var(g_rng);
          int flips_ref =  binom_noise_ref(g_rng);
          
          R.setVar(p, c, org_var - flips_var + flips_ref);
          R.setRef(p, c, coverage - R.getVar(p, c));
        }
        else
        {
          R.setVar(p, c, org_var);
          R.setRef(p, c, org_ref);
        }
      }
    }
  }
  
  return R;
}
