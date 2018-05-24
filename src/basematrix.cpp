/*
 * basematrix.cpp
 *
 *  Created on: 5-sep-2017
 *      Author: M. El-Kebir
 */

#include "basematrix.h"

BaseMatrix::BaseMatrix()
  : _m(0)
  , _k(0)
  , _n(0)
  , _indexToAnatomicalSite()
  , _anatomicalSiteToIndex()
  , _indexToSample()
  , _sampleToIndex()
  , _indexToCharacter()
  , _characterToIndex()
  , _sampleIndexToAnatomicalSiteIndex()
  , _anatomicalSiteIndexToSampleIndices()
{
}

BaseMatrix::BaseMatrix(const StringVector& indexToLocation,
                       const StringVector& indexToSample,
                       const StringVector& indexToCharacter,
                       const IntVector& sampleIndexToLocationIndex)
  : _m(indexToLocation.size())
  , _k(indexToSample.size())
  , _n(indexToCharacter.size())
  , _indexToAnatomicalSite(indexToLocation)
  , _anatomicalSiteToIndex()
  , _indexToSample(indexToSample)
  , _sampleToIndex()
  , _indexToCharacter(indexToCharacter)
  , _characterToIndex()
  , _sampleIndexToAnatomicalSiteIndex(sampleIndexToLocationIndex)
  , _anatomicalSiteIndexToSampleIndices(_m)
{
  for (int s = 0; s < _m; ++s)
  {
    const std::string& str = _indexToAnatomicalSite[s];
    _anatomicalSiteToIndex[str] = s;
  }
  
  for (int p = 0; p < _k; ++p)
  {
    const std::string& str = _indexToSample[p];
    _sampleToIndex[str] = p;
    int s = _sampleIndexToAnatomicalSiteIndex[p];
    _anatomicalSiteIndexToSampleIndices[s].insert(p);
  }
  
  for (int c = 0; c < _n; ++c)
  {
    const std::string& str = _indexToCharacter[c];
    _characterToIndex[str] = c;
  }
}
