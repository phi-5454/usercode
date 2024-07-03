#include "../interface/RpReco.h"

#include "../interface/Structures.h"
#include "../interface/RpPat.h"
#include "../interface/RpFit.h"
#include "../interface/RpEffic.h"
#include "../interface/RpPolygon.h"

#include "../interface/Helper.h"

#include <iomanip>

#define sqr(x) ((x)*(x))

using namespace std;

/*****************************************************************************/
RpReco::RpReco(int nCmTra, bool readRpAlign) : nCmTra(nCmTra)
{
  map<RpDet, map<vector<double>,int>> dummy;
  theRpFit   = new RpFit(dummy,nCmTra, 2); // read shifts, pattern fits
  theRpEffic = new RpEffic(2,nCmTra);      // read group effic

  //
  {
    int arm = 0;

    // near
    v_x[arm][0] = -2.20419654; L_x[arm][0] =   3104.19818816;
    v_y[arm][0] =  0.03239461; L_y[arm][0] = 238218.46794102;

    // far
    v_x[arm][1] = -1.88422776; L_x[arm][1] =   -522.47623562;
    v_y[arm][1] =  0.00750900; L_y[arm][1] = 271324.15847026;
  }

  //
  {
    int arm = 1;

    // near
    v_x[arm][0] = -2.24546659; L_x[arm][0] =    194.31957649;
    v_y[arm][0] =  0.01851263; L_y[arm][0] = 238326.01074279;

    // far
    v_x[arm][1] = -1.92329688; L_x[arm][1] =  -2950.79562493;
    v_y[arm][1] = -0.00829497; L_y[arm][1] = 271314.32190238;
  }

  // initialize corrections [arm][sta][rpt], legacy
  mx[0][1][0]=-0.465; mx[0][0][0]=-0.210;
  mx[1][0][0]= 0.167; mx[1][1][0]=-0.450;
  mx[0][1][1]=-0.081; mx[0][0][1]=-0.112;
  mx[1][0][1]= 0.373; mx[1][1][1]=-0.574;
  
  my[0][1][0]=-0.689; my[0][0][0]=-1.479;
  my[1][0][0]=-0.916; my[1][1][0]= 0.044;
  my[0][1][1]= 0.009; my[0][0][1]= 0.842;
  my[1][0][1]= 1.312; my[1][1][1]= 0.316;

  // additional shifts based on 0part TB and BT
  mx[0][0][0] -=  21e-3; // 1nT
  mx[0][0][1] -=  44e-3; // 1nB
  mx[0][1][0] -=   6e-3; // 1fT
  mx[0][1][1] -=  47e-3; // 1fB
  //
  mx[1][0][0] -=   8e-3; // 2nT
  mx[1][0][1] -=  47e-3; // 2nB
  mx[1][1][0] -=  -3e-3; // 2fT
  mx[1][1][1] -=  37e-3; // 2fB

  // additional shifts based on y* distribution
  for(int sta = 0; sta < 2; sta++)
  {
    int sign = (sta == 0 ? 1 : - 1);

    my[0][sta][0] -=  855e-3 / 70 * sign; // +- 12e-3 | 1xT
    my[0][sta][1] -=  750e-3 / 70 * sign; // +- 11e-3 | 1xB

    my[1][sta][0] -= 1540e-3 / 70 * sign; // +- 22e-3 | 2xT
    my[1][sta][1] -=  570e-3 / 70 * sign; // +-  8e-3 | 2xB
  }

  // read additional shifts from alignment
  if(readRpAlign)
  {
    ifstream file("../pars/rp_align.dat");
    cerr << Helper::col(2) << " reading roman pot alignment per run"
         << Helper::col();

    while(!file.eof())
    { 
      int run;
      file >> run;

      if(!file.eof())
      {
        cerr << ".";
        float d;

        for(int arm = 0; arm < 2; arm++)
        for(int sta = 0; sta < 2; sta++)
        for(int rpt = 0; rpt < 2; rpt++)
        { file >> d; dmx[arm][sta][rpt][run] = -d; }

        for(int arm = 0; arm < 2; arm++)
        for(int sta = 0; sta < 2; sta++)
        for(int rpt = 0; rpt < 2; rpt++)
        { file >> d; dmy[arm][sta][rpt][run] = -d; }

        for(int k = 0; k < 2; k++) // read uncertainties
        for(int arm = 0; arm < 2; arm++)
        for(int sta = 0; sta < 2; sta++)
        for(int rpt = 0; rpt < 2; rpt++)
        { file >> d; }
      }
    }
    cerr << " [done]" << endl;
  
    file.close();
  }
  else
    cerr << " not reading roman pot alignment per run" << endl;

  // histograms
  const string pre = "../out/rp/"+to_string(nCmTra)+"part";

  for(int arm = 0; arm < 2; arm++)
  for(int topo = 0; topo < nTopos; topo++)
  {
    const string post = to_string(arm)+"_"+topos[topo]+".his.gz";

    his_prot_loc[arm][topo].init(-1,1,100, -4,4,100,
        pre+"/prot/loc_"+post,false,true);
    his_prot_mom[arm][topo].init(-1,1,100, -1,1,100,
        pre+"/prot/mom_"+post,false,true);

    for(auto & run : runs)
      his_prot_loc_vs_run[arm][topo][run].init(-1,1,100, -4,4,100,
        pre+"/prot/loc_"+to_string(run)+"_"+post,false,true);
  }

  for(int topo = 0; topo < nTopos; topo++)
  {
    const string post = topos[topo]+".his";

    his_weight[topo].init(0,10,400, pre+"/weight_"+post);

    his_vtx_x[topo].init(-1,1,100, -1,1,100, pre+"/vtx/x_"+post);
    his_vtx_y[topo].init(-5,5,100, -5,5,100, pre+"/vtx/y_"+post);

    his_prot_mom_x[topo].init(-1,1,200, -1,1,200, pre+"/prot/mom_x_"+post);
    his_prot_mom_y[topo].init(-1,1,200, -1,1,200, pre+"/prot/mom_y_"+post);
  }

  for(int arm = 0; arm < 2; arm++)
  for(int sta = 0; sta < 2; sta++)
  for(int rpt = 0; rpt < 2; rpt++)
  for(int topo = 0; topo < nTopos; topo++)
  {
    const string post =
      to_string(arm)+to_string(sta)+to_string(rpt)+"_"+topos[topo]+".his.gz";

    his_hits_loc[arm][sta][rpt][topo].init(-2,2,200, -40,40,200,
        pre+"/hits/loc_"  +post,false,true);
    his_hits_loc_x[arm][sta][rpt][topo].init( -2, 2,800,
        pre+"/hits/loc_x_"+post,false,true);
    his_hits_loc_y[arm][sta][rpt][topo].init(-40,40,800,
        pre+"/hits/loc_y_"+post,false,true);

    for(auto & run : runs) 
    {
      his_hits_loc_vs_run[arm][sta][rpt][topo][run].init(-2,2,200, -40,40,200,
          pre+"/hits/loc_" +to_string(run)+"_"+post,false,true);

      his_hits_loc_x_vs_run[arm][sta][rpt][topo][run].init( -2, 2,800,
          pre+"/hits/loc_x_"+to_string(run)+"_"+post,false,true);

      his_hits_loc_y_vs_run[arm][sta][rpt][topo][run].init(-40,40,800,
          pre+"/hits/loc_y_"+to_string(run)+"_"+post,false,true);
    }
  }

  for(int arm = 0; arm < 2; arm++)
  for(int topo = 0; topo < nTopos; topo++)
  {
    const string post = to_string(arm)+"_"+topos[topo]+".his.gz";

    his_hits_nf_x[arm][topo].init( -2, 8,200,  -2, 8,200,
        pre+"/hits/nf_x_"+post,false,true);
    his_hits_nf_y[arm][topo].init(-40,40,200, -40,40,200,
        pre+"/hits/nf_y_"+post,false,true);

    for(auto & run : runs)
    {
      his_hits_nf_x_vs_run[arm][topo][run].init( -2, 8,200,  -2, 8,200,
          pre+"/hits/nf_x_"+to_string(run)+"_"+post,false,true);
      his_hits_nf_y_vs_run[arm][topo][run].init(-40,40,200, -40,40,200,
          pre+"/hits/nf_y_"+to_string(run)+"_"+post,false,true);

      his_hits_nf_dy_vs_run[arm][topo][run].init(-0.2,0.2,400,
         pre+"/hits/nf_dy_"+to_string(run)+"_"+post,false,true);
    }
  }

  //
  for(int a = 0; a < 2; a++)
  for(int s = 0; s < 2; s++)
  for(int r = 0; r < 2; r++)
  for(int u = 0; u < 2; u++)
  {
    RpDet det;
    det.arm = a; det.sta = s; det.rpt = r; det.uv  = u;

    char name[256];
    sprintf(name,"../out/rp/%dpart/occup/319311_%d%d%d_%d.his",
                 nCmTra, a,s,r, u);

    his_occup[det].init(0,nStrips,nStrips, -maxSlope,maxSlope,binSlope, name);
  }
}

/*****************************************************************************/
RpReco::~RpReco()
{
  ofstream fileTex("../an/rp_optics_"+to_string(nCmTra)+"part.tex");

  for(int arm = 0; arm < 2; arm++)
  for(int topo = 0; topo < nTopos; topo++)
  {
    double sig_xstar,sig_tstar, Ln,Lf;

    covariance(arm, his_hits_nf_x[arm][topo],
               sig_xstar,sig_tstar, Ln,Lf, topo);

    fileTex << fixed
      << " " << arm+1
      << " & " << topos[topo]
      << " & " << setprecision(0) << round(sig_xstar*1e+3)
      << " & " << setprecision(0) << round(sig_tstar*1e+6)
      << " & " << setprecision(0) << round(Ln/1e+1)*1e+1
      << " & " << setprecision(0) << round(Lf/1e+1)*1e+1
      << " & " << setprecision(arm==0?2:1)
                                  << round(Lf/Ln*1e+2)/1e+2
      << " \\\\ " << endl;
  }

  fileTex.close();
}

/*****************************************************************************/
void RpReco::covariance(int arm, Histo & his,
                        double & sig_xstar, double & sig_tstar,
                        double & Ln, double & Lf, int topo)
{
  const double & v1 = v_x[arm][0];
  const double & v2 = v_x[arm][1];

  double c11=0, c12=0, c22=0, vt=0, sum=0;
  const double det = 7000.; // mm

  for(int ix = 0; ix < his.axes[0].bins; ix++)
  for(int iy = 0; iy < his.axes[1].bins; iy++)
  {
    double x = his.axes[0].center(ix);
    double y = his.axes[1].center(iy);

    if(abs(x) < 0.7 && abs(y) < 0.7) // [mm] // PARAMETER
    {
      double val = his.get({ix,iy}).first;

      c11 += x*x * val;
      c12 += x*y * val;
      c22 += y*y * val;

      vt += sqr((v2*x - v1*y)/det) * val; // var(theta*)

      sum += val;
    }
  }

  c11 /= sum;
  c12 /= sum;
  c22 /= sum;
  vt  /= sum;

  // var(x*)
  double vx = (c11*c22 - sqr(c12))/(c22*sqr(v1) - 2*c12*v1*v2 + c11*sqr(v2));

  sig_xstar = sqrt(vx);
  sig_tstar = sqrt(vt);

  double L1_ = sqrt( (c11 - vx*sqr(v1)) / vt );
  double L2_ = sqrt( (c22 - vx*sqr(v2)) / vt );

  Lf = -L2_;
  Ln =  L1_;
}

/*****************************************************************************/
bool RpReco::getTopo(const vector<RpTrack> & rpTracks, int & topo)
{
  int nConf = 0;
  topo = -1;

  // put together code
if(0)
  {
    char code[4] = {'x','x','x','x'};

    for(auto & track : rpTracks)
    {
      RpDet det = track.det;
   
      int loc = -1;

      if(det.arm == 0 && det.sta == 0) loc = 1;
      if(det.arm == 0 && det.sta == 1) loc = 0;

      if(det.arm == 1 && det.sta == 0) loc = 2;
      if(det.arm == 1 && det.sta == 1) loc = 3;

      code[loc] = (det.rpt == 0 ? 'T' : 'B');
    }

mtx_cout.lock();
    cout << " " << code << endl;
mtx_cout.unlock();
  }

  // must have all four tracklets
  if(rpTracks.size() == 4)
  {
    bool v[2][2][2];
    for(short int arm = 0; arm < 2; arm++)
    for(short int sta = 0; sta < 2; sta++)
    for(short int rpt = 0; rpt < 2; rpt++)
      v[arm][sta][rpt] = false;
    
    for(auto & track : rpTracks)
    {
      RpDet det = track.det;

      const short int & arm = det.arm;
      const short int & sta = det.sta;
      const short int & rpt = det.rpt;
      
      v[arm][sta][rpt] = true;
    }

    // trigger
    bool lT = (v[0][1][0] && v[0][0][0]);
    bool lB = (v[0][1][1] && v[0][0][1]);

    bool rT = (v[1][0][0] && v[1][1][0]);
    bool rB = (v[1][0][1] && v[1][1][1]);
    
    if(lT && rB) { topo = 0; nConf++; }
    if(lB && rT) { topo = 1; nConf++; }
    if(lT && rT) { topo = 2; nConf++; }
    if(lB && rB) { topo = 3; nConf++; }
  }
  
  return (nConf == 1);
}

/*****************************************************************************/
void RpReco::localReco(int topo, int run, RpTrack & track,
                       map<RpDet,double> & locCy,
                       map<RpDet,double> & locCx)
{
  RpDet det = track.det;

  vector<Vector1> lpos(2);
  vector<double> eff(2); // u/v

  for(int uv = 0; uv < 2; uv++)
  {
    det.uv = uv;

    // copy
    vector<double> strip;
    for(int pla = 0; pla < nPlanes; pla++)
      strip.push_back(track.clus[uv][pla]);

    vector<double> orig = strip;

    RpPat theRpPat(4, -1); // don't read, nCmTra is not used
    vector<int> base_step = theRpPat.toPattern(strip); // strip mod!

    //
    LocalFit localFit;

    if(theRpFit->hasPattern(det,strip))
    { // already fitted
      const int & ibas = base_step[0];
      const int & base = base_step[1];
      const int & step = base_step[2];

      localFit = theRpFit->getLocalPosition(det,strip, ibas,base,step);
    }
    else
    { // something new (why?), fit now
      mtx_rpfit.lock();
      localFit = theRpFit->fitTracklet(det,orig);
      mtx_rpfit.unlock();
    }

    // local hit position
    lpos[uv] = {localFit.Cy, sqrt(localFit.Vy)}; 

    // get efficiency
//    eff[uv] = theRpEffic->getGroupEffic(run,det, localFit.Cy, localFit.Cx);
    eff[uv] = 1.; // FIXME

    locCy[det] = localFit.Cy;
    locCx[det] = localFit.Cx;
  } 

  track.weight = 1 / (eff[0] * eff[1]); // 1 / (tracklet uv[0] * tracklet uv[1])

  // get global hit position, overwrite!
  track.pos = theRpFit->getGlobalPosition(track.det, lpos);
}

/*****************************************************************************/
void RpReco::correctPos(int run, RpTrack & track)
{
  const auto & det = track.det;

  const auto & arm = det.arm;
  const auto & sta = det.sta;
  const auto & rpt = det.rpt;

  track.pos.x += mx[arm][sta][rpt] + dmx[arm][sta][rpt][run];
  track.pos.y += my[arm][sta][rpt] + dmy[arm][sta][rpt][run];
}

/*****************************************************************************/
void RpReco::getStar(int arm, const Vector2 pos[], RpRes & res)
{
  // 0 = near, 1 = far

  // determinant, 7000 mm
  double Dx = v_x[arm][0] * L_x[arm][1] -
              v_x[arm][1] * L_x[arm][0];

  double Dy = v_y[arm][0] * L_y[arm][1] -
              v_y[arm][1] * L_y[arm][0];

  // theta*
  {
    const double a1 =   v_x[arm][0]/Dx;
    const double a0 = - v_x[arm][1]/Dx;

    res.theta.x = a0*pos[0].x + a1*pos[1].x;

    const double b0 = 1/L_y[arm][0]/2;
    const double b1 = 1/L_y[arm][1]/2;

    res.theta.y = b0*pos[0].y + b1*pos[1].y;

    res.theta.cxx = sqr(a0)*pos[0].cxx + sqr(a1)*pos[1].cxx;
    res.theta.cyy = sqr(b0)*pos[0].cyy + sqr(b1)*pos[1].cyy;
    res.theta.cxy = a0* b0* pos[0].cxy + a1* b1* pos[1].cxy;
  }

  // x*, y*
  {
    const double a0 =   L_x[arm][1]/Dx;
    const double a1 = - L_x[arm][0]/Dx;

    res.pos.x = a0*pos[0].x + a1*pos[1].x;

    const double b0 =   L_y[arm][1]/Dy;
    const double b1 = - L_y[arm][0]/Dy;

    res.pos.y = b0*pos[0].y + b1*pos[1].y;

    res.pos.cxx = sqr(a0)*pos[0].cxx + sqr(a1)*pos[1].cxx;
    res.pos.cyy = sqr(b0)*pos[0].cyy + sqr(b1)*pos[1].cyy;
    res.pos.cxy = a0* b0* pos[0].cxy + a1* b1* pos[1].cxy;
  }
}

/*****************************************************************************/
void RpReco::globalReco(int topo, int run,
                        const vector<RpTrack> & rpTracks,
                              vector<PrTrack> & prTracks,
                        const double trkWeight,
                        double & rpWeight)
{
  //
  Vector2 pos[2][2]; // [arm][sta]
  rpWeight = 1.;

  // each tracklet
  for(auto & track : rpTracks)
  {
    RpDet det = track.det;

    const short int & arm = det.arm;
    const short int & sta = det.sta;

    pos[arm][sta] = track.pos;

    // Roman pots weight = Prod(proton weights) = Prod(tracklet weights)
    rpWeight *= track.weight;
  }

//  const double w = rpWeight * trkWeight; FIXME

  // each arm
  for(int arm = 0; arm < 2; arm++)
  {
    RpRes res;
    getStar(arm, pos[arm], res); // calculate

    // proton tracks
    PrTrack track;

       track.p.x = - BeamP      * res.theta.x; // p* - sign!
       track.p.y =   BeamP      * res.theta.y;

    // set rp track.pt here
    track.pt.x   = track.p.x; // copy
    track.pt.y   = track.p.y; // copy

    track.pt.cxx =   sqr(BeamP) * res.theta.cxx;
    track.pt.cxy = - sqr(BeamP) * res.theta.cxy;
    track.pt.cyy =   sqr(BeamP) * res.theta.cyy;

       track.p.z = (arm == 0 ? +BeamP : -BeamP);

       track.pos = res.pos; // x*, y*

    prTracks.push_back(track);
  }
}

/*****************************************************************************/
// reco and setp rpWeight
bool RpReco::reconstruct(Event & event)
{
  const int & topo = event.topo;
  const int & run  = event.run;

  vector<RpTrack> & rpTracks = event.rpTracks;
  vector<PrTrack> & prTracks = event.prTracks;

  // for his_occup
  map<RpDet,double> locCy,locCx;

  // re-reco and update rpTrack.pos
  for(auto & rpTrack : rpTracks)
  {
    localReco(topo,run, rpTrack, locCy,locCx);

    correctPos(run, rpTrack);
  }

  // global reco
  mtx_his.lock();
  globalReco(topo,run, rpTracks, prTracks, event.trkWeight,
                                           event.rpWeight); // get rpWeight
  mtx_his.unlock();

  // are protons within py limits?
  bool isWithin = false;

  if(fabs(prTracks[0].p.y) > minPy && fabs(prTracks[0].p.y) < maxPy &&
     fabs(prTracks[1].p.y) > minPy && fabs(prTracks[1].p.y) < maxPy)
    isWithin = true;

  return isWithin;
}

