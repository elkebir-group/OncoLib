/*
 * readmatrix.h
 *
 *  Created on: 5-sep-2017
 *      Author: M. El-Kebir
 */

#ifndef READMATRIX_H
#define READMATRIX_H

#include "utils.h"
#include "basematrix.h"
#include "frequencymatrix.h"

/// This class models a matrix of reads
class ReadMatrix : public BaseMatrix
{
public:
  /// Constructor
  ReadMatrix();
  
  /// Constructor
  ///
  /// @param indexToLocation Index to location label
  /// @param indexToSample Index to sample label
  /// @param indexToCharacter Index to character label
  /// @param sampleIndexToLocationIndex Sample index to location index
  ReadMatrix(const StringVector& indexToLocation,
             const StringVector& indexToSample,
             const StringVector& indexToCharacter,
             const IntVector& sampleIndexToLocationIndex);
  
  /// Return number of variant reads
  ///
  /// @param p Sample
  /// @param i Character
  int getVar(int p, int i) const
  {
    assert(0 <= p && p < _k);
    assert(0 <= i && i < _n);
    
    return _var[p][i];
  }
  
  /// Return number of reference reads
  ///
  /// @param p Sample
  /// @param i Character
  int getRef(int p, int i) const
  {
    assert(0 <= p && p < _k);
    assert(0 <= i && i < _n);
    
    return _ref[p][i];
  }
  
  /// Set number of variant reads
  ///
  /// @param p Sample
  /// @param i Character
  /// @param var Number of variant reads
  void setVar(int p, int i, int var)
  {
    assert(0 <= p && p < _k);
    assert(0 <= i && i < _n);
    
    _var[p][i] = var;
  }
  
  /// Set number of reference reads
  ///
  /// @param p Sample
  /// @param i Character
  /// @param ref Number of reference reads
  void setRef(int p, int i, int ref)
  {
    assert(0 <= p && p < _k);
    assert(0 <= i && i < _n);
    
    _ref[p][i] = ref;
  }
  
  /// Return frequency matrix using the specified alpha for computing
  /// the binomial proportion confidence intervals
  ///
  /// @param alpha Binomial proportion confidence interval
  /// @param threshold Variant read count threshold
  FrequencyMatrix toFrequencyMatrix(double alpha,
                                    int threshold) const;
  
  /// Return a read matrix whose entries are mutation clusters with pooled
  /// read counts
  ///
  /// @param clustering Cluster assignment of each SNV
  /// @param relabel Flag indicating whether clusters should receive new labels
  ReadMatrix poolReads(const IntMatrix& clustering,
                       bool relabel) const;
  
  /// Return a read matrix whose entries are down sampled from the current matrix
  ///
  /// @param nrSamplesPerAnatomicalSite Number of samples per anatomical site
  /// @param coverage Coverage
  /// @param purity Purity
  /// @param seqErrorRate Per base sequencing error rate
  /// @param fracSNVs Fraction of SNVs to retain
  ReadMatrix downSample(int nrSamplesPerAnatomicalSite,
                        int coverage,
                        double purity,
                        double seqErrorRate,
                        double fracSNVs) const;

private:
  /// Variant reads
  IntMatrix _var;
  /// Reference reads
  IntMatrix _ref;
  
  friend std::ostream& operator<<(std::ostream& out, const ReadMatrix& R);
  friend std::istream& operator>>(std::istream& in, ReadMatrix& R);
};

std::ostream& operator<<(std::ostream& out, const ReadMatrix& R);
std::istream& operator>>(std::istream& in, ReadMatrix& R);

#endif // READMATRIX_H
