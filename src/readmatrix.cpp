/*
 * readmatrix.cpp
 *
 *  Created on: 5-sep-2017
 *      Author: M. El-Kebir
 */

#include "readmatrix.h"
#include <boost/math/distributions/beta.hpp>

ReadMatrix::ReadMatrix()
  : BaseMatrix()
  , _var()
  , _ref()
{
}

ReadMatrix::ReadMatrix(const StringVector& indexToLocation,
                       const StringVector& indexToSample,
                       const StringVector& indexToCharacter,
                       const IntVector& sampleIndexToLocationIndex)
  : BaseMatrix(indexToLocation, indexToSample, indexToCharacter, sampleIndexToLocationIndex)
  , _var(IntMatrix(_k, IntVector(_n, 0)))
  , _ref(IntMatrix(_k, IntVector(_n, 0)))
{
}

ReadMatrix ReadMatrix::poolReads(const IntMatrix& clustering,
                                 bool relabel) const
{
  ReadMatrix newR(*this);
  newR._k = _k;
  newR._m = _m;
  newR._n = 0;
  newR._indexToCharacter.clear();
  newR._characterToIndex.clear();
  
  char buf[1024];
  
  // 1. infer character labeling
  IntVector oldCharacterToNewCharacter(_n, -1);
  int idx = 0;
  for (const IntVector& cluster : clustering)
  {
    std::string label;
    for (int i : cluster)
    {
      if (!label.empty())
        label += ";";
      
      label += indexToCharacter(i);
      oldCharacterToNewCharacter[i] = idx;
    }
    
    if (relabel)
    {
      snprintf(buf, 1024, "cluster_%d", idx + 1);
      label = buf;
    }
    
    newR._indexToCharacter.push_back(label);
    newR._characterToIndex[label] = idx;
    
    ++idx;
    ++newR._n;
  }
  
  // 2. pool reads
  newR._ref = IntMatrix(_k, IntVector(newR._n, 0));
  newR._var = IntMatrix(_k, IntVector(newR._n, 0));
  
  for (const IntVector& cluster : clustering)
  {
    for (int i : cluster)
    {
      for (int p = 0; p < _k; ++p)
      {
        int new_i = oldCharacterToNewCharacter[i];
        assert(new_i != -1);
        newR._ref[p][new_i] += _ref[p][i];
        newR._var[p][new_i] += _var[p][i];
      }
    }
  }
  
  return newR;
}

FrequencyMatrix ReadMatrix::toFrequencyMatrix(double alpha,
                                              int threshold) const
{
  FrequencyMatrix resF(*this);
  
  // compute confidence intervals
  for (int p = 0; p < _k; ++p)
  {
    for (int i = 0; i < _n; ++i)
    {
      int var = getVar(p, i);
      int ref = getRef(p, i);
      
      if (alpha < 0)
      {
        double vaf = (double)var / double(var + ref);
        resF.set(p, i, vaf, vaf);
      }
      else
      {
        boost::math::beta_distribution<> beta_dist(1 + var, 1 + ref);
        
        double f_lb = boost::math::quantile(beta_dist, alpha / 2);
        double f_ub = boost::math::quantile(beta_dist, 1 - alpha / 2);
        
        if (var < threshold || var <= (var + ref) * 0.01)
        {
          f_lb = f_ub = 0;
        }
        
  //      f_lb *= 2;
  //      f_ub *= 2;
        
        if (f_lb > 1)
        {
          f_lb = 1;
        }
        if (f_ub > 1)
        {
          f_ub = 1;
        }
        
        resF.set(p, i, f_lb, f_ub);
      }
    }
  }
  
  return resF;
}

ReadMatrix ReadMatrix::downSample(int nrSamplesPerAnatomicalSite,
                                  int coverage,
                                  double purity,
                                  double seqErrorRate,
                                  double fracSNVs) const
{
  // TODO: remove SNVs with 0 variant reads in all samples...
  ReadMatrix newR;
  newR._m = _m;
  newR._k = nrSamplesPerAnatomicalSite < 0 ? _k : nrSamplesPerAnatomicalSite * _m;
  newR._anatomicalSiteIndexToSampleIndices = IntSetVector(_m);
  
  int newP = 0;
  for (int s = 0; s < _m; ++s)
  {
    IntVector sampleIndices(anatomicalSiteIndexToSampleIndices(s).begin(),
                            anatomicalSiteIndexToSampleIndices(s).end());
    
    std::shuffle(sampleIndices.begin(), sampleIndices.end(), g_rng);
    
    int nrSamplesPerAnatomicalSite_s = nrSamplesPerAnatomicalSite < 0 ? anatomicalSiteIndexToSampleIndices(s).size() : nrSamplesPerAnatomicalSite;
    assert(nrSamplesPerAnatomicalSite_s <= sampleIndices.size());
    for (int pp = 0; pp < nrSamplesPerAnatomicalSite_s; ++pp)
    {
      int p = sampleIndices[pp];
      const std::string& pStr = _indexToSample[p];
      
      newR._sampleToIndex[pStr] = newP;
      newR._indexToSample.push_back(pStr);
      
      newR._anatomicalSiteIndexToSampleIndices[s].insert(newP);
      newR._sampleIndexToAnatomicalSiteIndex.push_back(s);
      
      ++newP;
    }
  }
  
  IntVector snvIndices;
  for (const std::pair<std::string, int> kv : _characterToIndex)
  {
    int i = kv.second;
    bool ok = false;
    for (const std::string& sample : _indexToSample)
    {
      if (newR._sampleToIndex.count(sample) == 1)
      {
        int p = _sampleToIndex.find(sample)->second;
        if (getVar(p, i) > 0)
        {
          ok = true;
          break;
        }
      }
    }
    
    if (ok)
    {
      snvIndices.push_back(i);
    }
  }
  
  const int new_n = std::min((int)(_n * fracSNVs), (int)snvIndices.size());
  std::shuffle(snvIndices.begin(), snvIndices.end(), g_rng);
  snvIndices.erase(snvIndices.begin() + new_n, snvIndices.end());
  newR._n = new_n;

  std::poisson_distribution<> poisson(coverage >= 0 ? coverage : 0);

  newR._var = IntMatrix(newR._k,
                        IntVector(new_n, 0));
  newR._ref = IntMatrix(newR._k,
                        IntVector(new_n, 0));
  
  newR._indexToAnatomicalSite = _indexToAnatomicalSite;
  newR._anatomicalSiteToIndex = _anatomicalSiteToIndex;
  newR._anatomicalSiteIndexToSampleIndices = IntSetVector(_m);
  newR._indexToCharacter = StringVector(new_n);
  for (int i = 0; i < new_n; ++i)
  {
    newR._indexToCharacter[i] = _indexToCharacter[snvIndices[i]];
    newR._characterToIndex[newR._indexToCharacter[i]] = i;
  }

  for (int p = 0; p < newR._k; ++p)
  {
    int pp = _sampleToIndex.find(newR.indexToSample(p))->second;
    for (int i = 0; i < new_n; ++i)
    {
      int old_i = snvIndices[i];
      if (coverage < 0)
      {
        newR._var[p][i] = _var[pp][old_i];
        newR._ref[p][i] = _ref[pp][old_i];
      }
      else
      {
        double vaf_pi = purity * double(_var[pp][old_i]) / double(_var[pp][old_i] + _ref[pp][old_i]);
        int newCoverage = poisson(g_rng);
        
        std::binomial_distribution<> binom(newCoverage, vaf_pi);
        int org_var = binom(g_rng);
        int org_ref = newCoverage - org_var;
        
        if (g_tol.nonZero(seqErrorRate))
        {
          std::binomial_distribution<> binom_noise_var(org_var,
                                                       seqErrorRate);
          std::binomial_distribution<> binom_noise_ref(org_ref,
                                                       seqErrorRate);
          
          int flips_var = binom_noise_var(g_rng);
          int flips_ref = binom_noise_ref(g_rng);
          
          newR._var[p][i] = org_var - flips_var + flips_ref;
          newR._ref[p][i] = newCoverage - newR._var[p][i];
        }
        else
        {
          newR._var[p][i] = org_var;
          newR._ref[p][i] = newCoverage - newR._var[p][i];
        }
      }
    }
  }
  
  return newR;
}

std::ostream& operator<<(std::ostream& out,
                         const ReadMatrix& R)
{
  out << R._m << " #anatomical sites" << std::endl;
  out << R._k << " #samples" << std::endl;
  out << R._n << " #characters" << std::endl;
  out << "#sample_index\tsample_label\tanatomical_site_index\t"\
         "anatomical_site_label\tcharacter_index\tcharacter_label\tref\tvar"
      << std::endl;
  for (int p = 0; p < R._k; ++p)
  {
    int s = R.sampleIndexToAnatomicalSiteIndex(p);
    const std::string sStr = R.indexToAnatomicalSite(s);
    
    const std::string& pStr = R.indexToSample(p);
    for (int c = 0; c < R._n; ++c)
    {
      const std::string& cStr = R.indexToCharacter(c);
      out << p << "\t" << pStr << "\t"
          << s << "\t" << sStr << "\t"
          << c << "\t" << cStr << "\t"
          << R._ref[p][c] << "\t" << R._var[p][c] << std::endl;
    }
  }
  
  return out;
}

std::istream& operator>>(std::istream& in, ReadMatrix& R)
{
  std::string line;
  getline(in, line);
  
  int m = -1;
  int k = -1;
  int n = -1;
  
  std::stringstream ss(line);
  ss >> m;
  
  if (m <= 0)
  {
    throw std::runtime_error("Error: m should be nonnegative");
  }
  
  getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> k;
  
  if (k < m)
  {
    throw std::runtime_error("Error: k should be at least m");
  }
  
  getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> n;
  
  if (n <= 0)
  {
    throw std::runtime_error("Error: n should be nonnegative");
  }
  
  R._m = m;
  R._k = k;
  R._n = n;
  
  R._indexToAnatomicalSite = StringVector(m);
  R._indexToSample = StringVector(k);
  R._indexToCharacter = StringVector(n);
  R._sampleIndexToAnatomicalSiteIndex = IntVector(k, -1);
  R._anatomicalSiteIndexToSampleIndices = IntSetVector(m);
  
  R._var = IntMatrix(k, IntVector(n, 0));
  R._ref = IntMatrix(k, IntVector(n, 0));
  
  std::vector<std::vector<bool> > present(k, std::vector<bool>(n, false));
  while (in.good())
  {
    getline(in, line);
    if (line == "" || line[0] == '#')
      continue;
    
    ss.clear();
    ss.str(line);
    
    int p = -1;
    std::string pStr;
    int s = -1;
    std::string sStr;
    int c = -1;
    std::string cStr;
    int ref = -1;
    int var = -1;
    
    ss >> p >> pStr >> s >> sStr >> c >> cStr >> ref >> var;
    
    if (!(0 <= p && p < k) || !(0 <= s && s < m) || !(0 <= c && c < n)
        || ref < 0 || var < 0)
    {
      throw std::runtime_error("Invalid character ("
                               + boost::lexical_cast<std::string>(p)
                               + ","
                               + boost::lexical_cast<std::string>(c)
                               + ")");
    }
    
    if (present[p][c])
    {
      throw std::runtime_error("Duplicate character ("
                               + boost::lexical_cast<std::string>(p)
                               + ","
                               + boost::lexical_cast<std::string>(c)
                               + ")");
    }
    
    R._indexToAnatomicalSite[s] = sStr;
    R._anatomicalSiteToIndex[sStr] = s;
    R._indexToSample[p] = pStr;
    R._sampleToIndex[pStr] = p;
    R._indexToCharacter[c] = cStr;
    R._characterToIndex[cStr] = c;
    R._var[p][c] = var;
    R._ref[p][c] = ref;
    R._sampleIndexToAnatomicalSiteIndex[p] = s;
    R._anatomicalSiteIndexToSampleIndices[s].insert(p);
    
    present[p][c] = true;
  }
  
  for (int c = 0; c < R._n; ++c)
  {
    const std::string& str = R._indexToCharacter[c];
   
    StringVector s;
    boost::split(s, str, boost::is_any_of(";"));
    R._characterMultiplicity.push_back(s.size());
  }
  
  return in;
}
