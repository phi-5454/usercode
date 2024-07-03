#ifndef _Parameters_h_
#define _Parameters_h_

#include <vector>
#include <cmath>
#include <string>

///////////////////////////////////////////////////////////////////////////////
// beam
const double BeamP = 6500; // GeV

///////////////////////////////////////////////////////////////////////////////
// runs
const std::vector<int> runs =
    { 319104, 319124, 319125, 319159, 319160,
      319174, 319175, 319176, 319177, 319190,
      319222, 319223, 319254, 319255, 319256,
      319260, 319262, 319263, 319264, 319265,
      319266, 319267, 319268, 319270, 319300, 319311 } ; // 26 runs

///////////////////////////////////////////////////////////////////////////////
// roman pots
const double maxSlope = 0.4; // PARAMETER
const int    binSlope = 40;  // PARAMETER

//
enum                                    { TB=0, BT=1, TT=2, BB=3, nTopos};
const std::vector<std::string> topos = { "TB", "BT", "TT", "BB" };

//
const double maxDphi = M_PI, minKt = 0.20, maxKt = 0.80; // PARAMETER
const int      nDphi = 18  ,   nKt = 12;                 // PARAMETER
const int                      bKt = 6; // for pwav, spherical harmonics

// min/max p_y of scattered protons
const double minPy = 0.175; // GeV
const double maxPy = 0.670; // GeV

const int Empty = -99; // empty 

// max(hat(t),hat(u))
const double minSha = -3, maxSha = 1;
const int      nSha = 40;

// particles
enum                                   { pion, kaon, prot, nParts, unknown };
const std::vector<std::string> parts = { "pi", "ka", "pr" };
const std::vector<double> mass =       { 0.139570, 0.493677, 0.938272, 0, 0 };

// tracks: eta [-3:3, bin 0.1], pt [0:2, bin 50 MeV]
const double maxEta = 3., maxPt = 2., maxPsi = M_PI; // PARAMETER
const int      nEta = 60,   nPt = 40,   nPsi = 36;   // PARAMETER

// event selection
enum { nothing, elastic, signal, sideband };

// track pair
// accepted y region
const std::vector<double> maxResoY = {2.0, 2.0, 2.0}; // PARAMETER

// mass [0-4, bin 20 MeV], (cosTheta,phi) [bin 10x10]
const double maxMass = 4., maxCosTheta = 1., maxPhi = M_PI; // PARAMETER
const int      nMass = 200,  nCosTheta = 10,   nPhi = 10;   // PARAMETER

const double minTkEff = 0.1; // minimal one-track tracking effic PARAMETER

const int nPidMcEvents = 1;

///////////////////////////////////////////////////////////////////////////////
// pwav
const std::vector<std::pair<int,int>> lms =
  {
    std::pair<int,int>(0,0),
    //
    std::pair<int,int>(1,0),
    std::pair<int,int>(1,1),
    //
    std::pair<int,int>(2,0),
    std::pair<int,int>(2,1),
    std::pair<int,int>(2,2),
    //
    std::pair<int,int>(3,0),
    std::pair<int,int>(3,1),
    std::pair<int,int>(3,2),
    std::pair<int,int>(3,3),
    //
    std::pair<int,int>(4,0),
    std::pair<int,int>(4,1),
    std::pair<int,int>(4,2),
    std::pair<int,int>(4,3),
    std::pair<int,int>(4,4),
  };

const int nLms = 1+2+3+4+5; // = lms.size()

enum { do_separate, go_up, go_down };

#endif
