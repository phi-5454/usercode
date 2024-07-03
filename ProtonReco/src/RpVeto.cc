#include "../interface/RpVeto.h"

#include "../interface/Structures.h"
#include "../interface/Helper.h"

#include <iostream>
#include <cmath>
#include <map>

using namespace std;

const int nEvents_rp_eff = 4e+8; // PARAMETER
const double minRpEff = 1e-2;    // PARAMETER

/*****************************************************************************/
// flag =-1 : do nothing
// flag = 1 : collect for veto
// flag = 0 : calculate angular efficiency
// flag = 2 : read and use roman pot efficiency (accep and !veto)
RpVeto::RpVeto(int flag, int nCmTra) : flag(flag)
{
  string base = "../out/rp/"+to_string(nCmTra)+"part/";

  if(flag == 1) // collect for veto
  {
    his_veto_effic.init(-1,1,400, -1,1,400, base+"veto_py.his");
    den_veto_effic.init(-1,1,400, -1,1,400, "");

    for(int topo = 0; topo < nTopos; topo++)
    {
      his_py12[topo].init(0.1,0.5,200,
                          0.1,0.5,200,
                          base+"py12_"+topos[topo]+".his");

      his_py12_pass[topo].init(0.1,0.5,200,
                               0.1,0.5,200,
                               base+"py12_"+topos[topo]+"_pass.his");

      his_py12_veto[topo].init(0.1,0.5,200,
                               0.1,0.5,200,
                               base+"py12_"+topos[topo]+"_veto.his");
    }
  }

  if(flag == 0) // read, calculate angular efficiency
  {
    cerr << Helper::col(2) << " reading roman pots !veto efficiency"
         << Helper::col();

    his_veto_effic.init(-1,1,400, -1,1,400, base+"veto_py.his");
    his_veto_effic.read();
 
    cerr << " [done]" << endl;
  }

  if(flag == 2) // read angular efficiency 
  {
    const double delKt = (maxKt-minKt)/3;

    cerr << Helper::col(2) << " reading roman pot angular efficiency"
         << Helper::col();
 
    for(int topo = 0; topo < nTopos; topo++)
    {
      cerr << ".";
      Histo & h = his_eff[topo];
      h.init(0,maxDphi,nDphi, minKt-delKt,maxKt+delKt,(5*nKt)/3,
                              minKt-delKt,maxKt+delKt,(5*nKt)/3,
             "../out/rp/eff_"+topos[topo]+".his");

      h.read();
    }
    cerr << " [done]" << endl;
  }
}

/*****************************************************************************/
RpVeto::~RpVeto()
{
  if(flag == 1) // normalize
    his_veto_effic.div(den_veto_effic);
}

/*****************************************************************************/
short int RpVeto::collectRegions(const vector<short int> & fixed)
{
  // collect regions
  map<short int,int> regs;
  for(short int pla = 0; pla < nPlanes; pla++)
  if(fixed[pla] != -99)
  {
    short int reg = fixed[pla] / wGroup;
    regs[reg]++;
  }

  // find region with highest occupancy
  short int mreg = -1;
  int rmax = -1;

  for(auto & r : regs)
  if(r.second > rmax)
  {
    mreg = r.first;
    rmax = r.second;
  }

  // region with max counts
  return mreg;
}

/*****************************************************************************/
// for re-applying on data
// for calculating efficiency tables
bool RpVeto::elasticVeto(const vector<RpTrack> & rpTracks, int sta,
                         bool print) // [4]
{
  short int regs[2][2]; // [topo][uv]

  for(auto & track : rpTracks)
  if(track.det.sta == sta) // far
  for(short int uv = 0; uv < 2; uv++)
  {
    vector<short int> strip(nPlanes);
    for(short int pla = 0; pla < nPlanes; pla++)
      strip[pla] = int(track.clus[uv][pla]); // rounding?

    // get regions
    const short int & arm = track.det.arm;
    regs[arm][uv] = collectRegions(strip);
  }

  // t-map
  bool veto = true;
  double y[2];

  const vector<int> c = {-1,1}; // PARAMETERS for {TB,BT}
  int topo = (rpTracks[0].det.rpt == 0 ? TB : BT);

  for(short int arm = 0; arm < 2; arm++)
  {
    veto &= (abs(regs[arm][0] - regs[arm][1] - c[topo]) <= 1 && // PARAMETER
                 regs[arm][0] + regs[arm][1] >= 18);

    y[arm] = 25 - (regs[arm][0] + regs[arm][1]);


    if(print)
      cout << " " << regs[arm][0] - regs[arm][1] - c[topo]
           << " " << regs[arm][0] + regs[arm][1]
           << " " << y[arm];
  }

  veto &= (y[0] - y[1] == 0);

  if(print)
    cout << " " << (veto ? "veto" : "ok" );

  return veto;
}

/*****************************************************************************/
// same for both arms; for luminosity check
bool RpVeto::elasticMask(const vector<PrTrack> & prTracks, int topo)
{
  bool mask = false;

  if(topo == TB || topo == BT)
  {
    const double w = 0.070;
    const double d = 0.010;

    const double & p1y = prTracks[0].p.y;
    const double & p2y = prTracks[1].p.y;

    for(int i = 2; i <= 5; i++)
    {
      if(topo == TB)
      if(p1y >=  0.015-d + i*w && p1y <   0.015+d + (i+1)*w)
      if(p2y <  -0.025+d - i*w && p2y >= -0.025-d - (i+1)*w)
      { mask = true; break; }

      if(topo == BT)
      if(p1y <  -0.025+d - i*w && p1y >= -0.025-d - (i+1)*w)
      if(p2y >=  0.015-d + i*w && p2y <   0.015+d + (i+1)*w)
      { mask = true; break; }
    }
  }

  return mask;
}

/*****************************************************************************/
bool RpVeto::isAccepted(double py)
{
  return (fabs(py) > minPy &&
          fabs(py) < maxPy);
}

/*****************************************************************************/
int RpVeto::getTopology(double p1y, double p2y)
{
  if(p1y >= 0 && p2y <  0) return TB;
  if(p1y <  0 && p2y >= 0) return BT;
  if(p1y >= 0 && p2y >= 0) return TT;
  if(p1y <  0 && p2y <  0) return BB;

  exit(1);
}

/*****************************************************************************/
double RpVeto::getAngularCoverage(int topo, const vector<float> & dphip1tp2t)
{
  return his_eff[topo].val(dphip1tp2t);
}


