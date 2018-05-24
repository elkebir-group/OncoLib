/*
 * mutclonetree.h
 *
 *  Created on: 21-may-2018
 *      Author: M. El-Kebir
 */

#ifndef MUTCLONETREE_H
#define MUTCLONETREE_H

#include "utils.h"
#include "clonetree.h"
#include "frequencymatrix.h"
#include "readmatrix.h"

/// This class models a clone tree with mutations
class MutCloneTree : public CloneTree
{
public:
  /// Default constructor
  MutCloneTree();
  
  /// Copy constructor
  MutCloneTree(const MutCloneTree& other);
  
  /// Constructor
  ///
  /// @param T Directed graph
  /// @param root Root node
  /// @param id Node identifier
  /// @param l Leaf label
  MutCloneTree(const Digraph& T,
               Node root,
               const StringNodeMap& id,
               const StringNodeMap& l);
  
  /// Assignment operator
  MutCloneTree& operator =(const MutCloneTree& other);
  
  /// Read leaf labeling
  ///
  /// @param in Input stream
  bool readLeafLabeling(std::istream& in);
  
  /// Write leaf labeling
  ///
  /// @param out Output stream
  void writeLeafLabeling(std::ostream& out) const;
  
  /// Return mixture proportion
  ///
  /// @param u Leaf
  const DoubleVector& getMixtureProportions(Node u) const
  {
    return _U[u];
  }
  
  /// Set mixture proportion
  ///
  /// @param u Leaf
  /// @param proportions Mixture proportions
  void setMixtureProportions(Node u, const DoubleVector& proportions)
  {
    _U[u] = proportions;
  }
  
  /// Return mutation state
  ///
  /// @param u Leaf
  /// @param c Character
  bool hasMutation(Node u, const std::string& c) const
  {
    return _b[u].count(c) == 1;
  }
  
  /// Return mutations
  ///
  /// @param u Leaf
  const StringSet& getMutations(Node u) const
  {
    return _b[u];
  }
  
  /// Return mutations node map
  const StringSetNodeMap& getMutationsNodeMap() const
  {
    return _b;
  }
  
  /// Set mutations
  ///
  /// @param u Leaf
  /// @param b Mutation set
  void setMutations(Node u, const StringSet& b)
  {
    _b[u] = b;
  }
  
  /// Add mutation
  ///
  /// @param u Leaf
  /// @param c Character
  void addMutation(Node u, const std::string& c)
  {
    _b[u].insert(c);
  }
  
  /// Remove mutation
  ///
  /// @param u Leaf
  /// @param c Character
  void removeMutation(Node u, const std::string& c)
  {
    _b[u].erase(c);
  }
  
  /// Print tree in DOT format using the given color map
  ///
  /// @param out Output stream
  /// @param colorMap Color map
  virtual void writeDOT(std::ostream& out,
                        const StringToIntMap& colorMap) const;
  
  /// Print tree in DOT format using the given vertex labeling and color map
  ///
  /// @param out Output stream
  /// @param lPlus Labeling of each node by a location
  /// @param colorMap Color map
  virtual void writeDOT(std::ostream& out,
                        const StringNodeMap& lPlus,
                        const StringToIntMap& colorMap) const;
  
  /// Generate frequency matrix
  FrequencyMatrix getFrequencies() const;
  
  /// Generate read matrix
  ///
  /// @param sequencingDepth Depth of sequencing
  /// @param sequencingErrorRate Nucleotide substition error rate
  /// @param purity Purity
  /// @param ploidy Ploidy
  ReadMatrix getReads(double purity,
                      int sequencingDepth,
                      double sequencingErrorRate,
                      int ploidy) const;
  
protected:
  /// Mixture proportion in location L(T) -> [0,1]^k
  DoubleVectorNodeMap _U;
  /// Mutation vector L(T) -> {0,1}^n
  StringSetNodeMap _b;
};

#endif // MUTCLONETREE_H
