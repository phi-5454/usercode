#ifndef _Commons_h_
#define _Commons_h_

#include <vector>
#include <sstream> // for ostringstream
#include <tuple>   // for tie

#include "../interface/Vectors.h"

//
// roman pots
const short int nPlanes =   5;
const short int nStrips = 512;
const short int wGroup  =  32;

//
struct RpDet {
  // arm (left/right); sta (near/far); rpt (top/bottom); u/v; plane (0-4);
  // strip
  short int arm, sta, rpt, uv=-1, pla=-1, str=-1;

  std::vector<double> clus; // for printing only

  std::string print() const {
    std::ostringstream os;

    os << arm<<"|"<<sta<<"|"<<rpt;

    if(uv !=-1) os <<"|u"<<uv ;
    if(pla!=-1) os <<"|p"<<pla;
    if(str!=-1) os <<"|r"<<str;

    return os.str();
  }

  bool operator <(const RpDet & y) const {
  return std::tie(  arm,   sta,   rpt,   uv,   pla,   str) <
         std::tie(y.arm, y.sta, y.rpt, y.uv, y.pla, y.str);
  }
};

//
struct RpRes {
  Vector2 pos;   // star
  Vector2 theta; // star
};

//
struct RpTrack { // tracklet in a roman pot, has u and v components
  RpDet det;
  Vector2 pos; // local position

  double clus[2][nPlanes]; // [uv][pla]

  double weight;
};

//
struct PrTrack // at IP
{
  Vector3 p;
  Vector2 pt;
  Vector2 pos;
};

//
struct CmVertex
{
  Vector1 x,y,z;

  double chi2;
  int ndf; 
};

// CmTrack
struct CmTrack
{
  Vector3 p;
  Vector2 pt;
  int q;

  double eps,sigma;

  Vector1 dt,dz;
  double chi2;
  int ndf;

  double sigpt;
};

//
struct Event
{
  int run, ls, bx;

  std::vector<RpDet>   rpDets;
  std::vector<RpTrack> rpTracks;

  std::vector<PrTrack> prTracks;

  CmVertex cmVertex;
  std::vector<CmTrack> cmTracks;

  int sign;
  double rpWeight;
  double trkWeight;

  int cat;  // category: nothing, elastic, signal, sideband
  int type; // central hadron: unknown, pion, kaon, prot
  int topo; // topology: TB, BT, TT, BB
};

#endif
