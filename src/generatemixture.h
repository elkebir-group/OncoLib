/*
 *  generatemixture.h
 *
 *   Created on: 22-may-2018
 *       Author: M. El-Kebir
 */

#ifndef GENERATEMIXTURE_H
#define GENERATEMIXTURE_H

#include "mutclonetree.h"
#include "utils.h"
#include "frequencymatrix.h"

class GenerateMixture
{
public:
  /// Constructor
  ///
  /// @param T Clone tree
  GenerateMixture(const MutCloneTree& T);
  
  /// Generate mixture
  ///
  /// @param k Number of samples per location
  /// @param partition Indicates whether clones per location should simply be paritioned.
  MutCloneTree generate(int k,
                        bool partition) const;
  
private:
  typedef std::map<std::string, NodeSet> StringToNodeSetMap;
  typedef std::map<int, NodeSet> IntToNodeSetMap;
  typedef std::map<int, int> IntToIntMap;
  
  MutCloneTree constructCloneTree(const NodeSet& sampledLeaves,
                                  const DoubleVectorNodeMap& sampleProportions) const;
  
private:
  /// Clone tree
  const MutCloneTree& _T;
};

#endif // GENERATEMIXTURE_H
