/*
 * simulatemain.cpp
 *
 *  Created on: 25-aug-2017
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "simulation.h"
#include <lemon/arg_parser.h>
#include <fstream>

int main(int argc, char** argv)
{
  int seed = 0;
  double K = 5e4;
  double mutFreqThreshold = 0.05;
  double migrationRate = 1e-6;
  double mutationRate = 0.1;
  double driverProb = 2e-7;
  int maxNrAnatomicalSites = 3;
  std::string filenameColorMap;
  std::string outputDirectory;
  int pattern = 0;
  int N = -1;
  int coverage = 200;
  int nrSamplesPerAnatomicalSite = 2;
  int nrSamplesPrimary = 2;
  double seqErrorRate = 0;
  double purity = 1;
  bool verbose = false;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("s", "Random number generator seed (default: 0)", seed)
    .refOption("c", "Color map file", filenameColorMap)
    .refOption("C", "Target coverage (default: 200)", coverage)
    .refOption("D", "Driver probability (default: 1e-7)", driverProb)
    .refOption("E", "Per base sequencing error rate (default: 0)", seqErrorRate)
    .refOption("P", "Purity (default: 1)", purity)
    .refOption("K", "Carrying capacity (default: 5e4)", K)
    .refOption("k", "Number of samples per anatomical site (default: 2)", nrSamplesPerAnatomicalSite)
    .refOption("kP", "Number of samples for the primary tumor (default: 2)", nrSamplesPrimary)
    .refOption("f", "Mutation frequency threshold (default: 0.05)", mutFreqThreshold)
    .refOption("mig", "Migration rate (default: 1e-6)", migrationRate)
    .refOption("mut", "Mutation rate (default: 0.1)", mutationRate)
    .refOption("N", "Number of successful simulations (default: -1)", N)
    .refOption("m", "Maximum number of detectable anatomical sites (default: 3)", maxNrAnatomicalSites)
    .refOption("p", "Allowed migration patterns:\n"\
               "       0 : mS (default)\n"\
               "       1 : mS, S\n" \
               "       2 : mS, S, M\n" \
               "       3 : mS, S, M, R", pattern)
    .refOption("v", "Verbose", verbose)
    .refOption("o", "Output directory", outputDirectory);
  ap.parse();
  
  if (Simulation::run(outputDirectory, filenameColorMap,
                      seed, K, mutFreqThreshold,
                      migrationRate, mutationRate, driverProb,
                      maxNrAnatomicalSites, pattern, N,
                      coverage, nrSamplesPerAnatomicalSite, nrSamplesPrimary,
                      seqErrorRate, purity, verbose) >= 1)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}
