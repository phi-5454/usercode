#ifndef _RpVeto_h_
#define _RpVeto_h_

#include <vector>

#include "../interface/Parameters.h"
#include "../interface/Histo.h"

struct RpTrack;
struct PrTrack;

//
class RpVeto
{
 public:
  RpVeto(int flag, int nCmTra);
  ~RpVeto();

  bool elasticVeto(const std::vector<RpTrack> & rpTracks, int sta, bool print);

  void combineTracks(const std::vector<RpTrack> & rpTracks,
                     const std::vector<PrTrack> & prTracks,
                     int topo, double rpWeight);

  void checkVeto(const std::vector<RpTrack> & rpTracks,
                 const std::vector<PrTrack> & prTracks,
                 int topo, double rpWeight);

  bool   elasticMask(const std::vector<PrTrack> & prTracks, int topo);

  void calculateAngularEfficiency();
  double      getAngularCoverage(int topo,
                                 const std::vector<float> & dphip1tp2t);

 private:
  short int collectRegions(const std::vector<short int> & fixed);

  bool isAccepted(double py);
  int getTopology(double p1y, double p2y);

  // for mixed events (elastic veto)
  std::vector<RpTrack> rpStack[2]; // [top or bottom]
  std::vector<PrTrack> prStack[2]; // [top or bottom]

  Histo his_veto_effic, den_veto_effic;

  Histo his_py12[nTopos], his_py12_veto[nTopos], his_py12_pass[nTopos];

  Histo his_eff[nTopos];

  int flag;

   std::mutex mtx_veto, mtx_check;
};

#endif
