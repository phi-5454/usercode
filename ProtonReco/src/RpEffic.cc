#include "../interface/RpEffic.h"

#include "../interface/Parameters.h"

#include "../interface/Helper.h"
#include "../interface/gzstream.h"
#include "../interface/RpPat.h"
#include "../interface/RpFit.h"
#include "../interface/Histo.h"

using namespace std;

/*****************************************************************************/
// flag = 0 : collect and write
// flag = 1 : read strip efficiency
// flag = 2 : reag group efficiency
RpEffic::RpEffic(const int flag, const int nCmTra) :
                      flag(flag), nCmTra(nCmTra)
{
  // collect
  if(flag == 0)
  {
    theRpPat = new RpPat(1, nCmTra);

    const map<RpDet, map<vector<double>,int>> dummy;
    theRpFit = new RpFit(dummy, nCmTra,2);
  }

  // use
  if(flag == 1)
    readStripEffic();

  //
/*
  if(flag == 2)
    readGroupEffic(); // FIXME
*/
}

/*****************************************************************************/
RpEffic::~RpEffic()
{
  // collect
  if(flag == 0)
  {
    writeStripEffic();

    delete theRpPat;
  }
}

/*****************************************************************************/
bool RpEffic::collectStrip(const vector<RpTrack> & rpTracks, const int & run)
{
  // take each
  for(auto & track : rpTracks)
  {
    RpDet det = track.det;

    // both orientation
    for(short int uv = 0; uv < 2; uv++)
    {
      det.uv = uv;

      // copy
      vector<double> str;
      for(short int pla = 0; pla < nPlanes; pla++)
        str.push_back(track.clus[uv][pla]);

      theRpPat->collectForEffic(det,str, theRpFit, nFound[run],nEmpty[run]);
    }
  }

  return true;
}

/*****************************************************************************/
void RpEffic::writeStripEffic(const int & run)
{
  char name[256];
  sprintf(name, "../out/rp/%dpart/stripEffic/%d_%d.dat.gz", nCmTra, run,nCmTra);
  ogzstream fileEff(name);

  //
  RpDet det;
  for(det.arm = 0; det.arm < 2; det.arm++)
  for(det.sta = 0; det.sta < 2; det.sta++)
  for(det.rpt = 0; det.rpt < 2; det.rpt++)
  for(det.uv  = 0; det.uv  < 2; det.uv++ )
  for(det.pla = 0; det.pla < nPlanes; det.pla++)
  for(det.str = 0; det.str < nStrips; det.str++)
  {
    const string key = det.print();

    const double & a = nFound[run][det];
    const double & b = nEmpty[run][det];

    double val = a/(a+b);
    double sig = sqrt(a*b/pow(a+b,3.)) + 0.5/(a+b);

    if(a == 0) val = 0;

    fileEff << " " << det.print() << " " << val << " " << sig << endl;
  }

  fileEff.close();
}

/*****************************************************************************/
void RpEffic::writeStripEffic()
{
  cerr << Helper::col(5) << " writing roman pots strip efficiency"
       << Helper::col();

  for(auto & run : runs)
  {
    writeStripEffic(run);
    cerr << ".";
  }

  cerr << " [done]" << endl;
}

/*****************************************************************************/
void RpEffic::readStripEffic(const int & run)
{
  char name[256];
  sprintf(name,"../out/rp/%dpart/stripEffic/%d_%d.dat.gz",nCmTra, run,nCmTra);
  igzstream fileEff(name);

  string line;

  while(getline(fileEff,line))
  {
    const vector<string> arg = Helper::parse(line," ");

    string reg = arg[0];

    double val = stof(arg[1]);
    double sig = stof(arg[2]);

    // (0 -nan -> 1) did not have enough hits
    if(val == 0 && std::isnan(sig)) val = 1;

    stripEffic[run][reg] = val;
  }

  fileEff.close();
}

/*****************************************************************************/
void RpEffic::readStripEffic()
{
  cerr << Helper::col(2)
       << " reading roman pots strip efficiencies"
       << Helper::col();

  for(auto & run : runs)
  {
    readStripEffic(run);
    cerr << ".";
  }

  cerr << " [done]" << endl;
}

/*****************************************************************************/
double RpEffic::getStripEffic(int run, const RpDet & det)
{
  return stripEffic[run][det.print()];
}

/*****************************************************************************/
void RpEffic::readGroupEffic(const int & run)
{
  for(int a = 0; a < 2; a++)
  for(int s = 0; s < 2; s++)
  for(int r = 0; r < 2; r++)
  for(int u = 0; u < 2; u++)
  {
    RpDet det;
    det.arm = a;
    det.sta = s;
    det.rpt = r;
    det.uv  = u;

    char name[256];
    sprintf(name,"../out/rp/%dpart/groupEffic/%d_%d%d%d_%d.his.gz",
                 2, run,a,s,r, u); // nCmTra = 2 (!) FIXME

    his_groupEffic[run][det].init(0,nStrips,nStrips,
                                  -maxSlope,maxSlope,binSlope, name,
                                  false, true); // zip
    his_groupEffic[run][det].read();
  }
}

/*****************************************************************************/
void RpEffic::readGroupEffic()
{
  cerr << Helper::col(2)
       << " reading roman pots group efficiencies (2 tracks)"
       << Helper::col();

  for(auto & run : runs)
  {
    readGroupEffic(run);
    cerr << ".";
  }

  cerr << " [done]" << endl;
}

/*****************************************************************************/
double RpEffic::getGroupEffic(int run, const RpDet & det,
                            double Cy, double Cx)
{
  double eff  = his_groupEffic[run][det].val({float(Cy),float(Cx)});

  // if it was outside (0,nStrips) or (-maxSlope,maxSlope)
  if(eff == 0) eff = 1;

  return eff;
}

