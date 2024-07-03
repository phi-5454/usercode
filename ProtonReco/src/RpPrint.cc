#include "../interface/RpPrint.h"

#include "../interface/Structures.h"

#include <iostream>

using namespace std;

enum { black=0, red=31, green=32, yellow=33, blue=34 };

/*****************************************************************************/
string RpPrint::col(int col = black) // 2print
{
  ostringstream os; os << "\033[" << col << "m"; return os.str();
}

/*****************************************************************************/
void RpPrint::print(bool hasDet, bool hasTra, int uv) // 2print
{
  if(!hasDet && !hasTra) cerr << " ";

  if( hasDet && !hasTra) cerr << (uv == 0 ? col(yellow) +"*" : col(green)+"*");
  if( hasDet &&  hasTra) cerr << (uv == 0 ? col(red)    +"x" : col(blue) +"+");

  if(!hasDet &&  hasTra) exit(1);

  cerr << col();
}

/*****************************************************************************/
void RpPrint::print(const vector<RpDet>   & rpDets,
                    const vector<RpTrack> & rpTracks) // 2print
{
  cerr << endl;
  const int nRegions = nStrips / wGroup;

  // rpTracks -> hit
  int hitDet[2][2][2][2*nPlanes][nRegions];
  int hitTra[2][2][2][2*nPlanes][nRegions];

  for(int arm = 0; arm < 2; arm++)
  for(int sta = 0; sta < 2; sta++)
  for(int rpt = 0; rpt < 2; rpt++)
  for(int pla = 0; pla < 2*nPlanes ; pla++)
  for(int reg = 0; reg < nRegions; reg++)
  {
    hitDet[arm][sta][rpt][pla][reg] = 0;
    hitTra[arm][sta][rpt][pla][reg] = 0;
  }

  for(auto & det : rpDets)
  for(auto & strip : det.clus)
  {
    int reg = int(strip) / wGroup;
    hitDet[det.arm][det.sta][det.rpt][det.uv + 2*det.pla][reg]++;
  }

  for(auto & track : rpTracks)
  {
    const RpDet & det = track.det;

    for(int uv = 0; uv < 2; uv++)
    for(int pla = 0; pla < nPlanes; pla++)
    {
      int istrip = int(track.clus[uv][pla]);

      if(istrip != -99)
      {
        int reg = istrip / wGroup;
        hitTra[det.arm][det.sta][det.rpt][uv + 2*pla][reg]++;
      }
    }
  }

  // print
  string dash(27,'-'), dots(27,'=');

  vector<pair<int,int>> pots = { {0,1}, {0,0}, {1,0}, {1,1} };

  cerr << " " << dash << "    " << dash << endl;
  for(int rpt = 0; rpt < 2; rpt++)
  {
    for(int reg = 0; reg < nRegions; reg++)
    {
      cerr << " | ";

      for(auto & pot : pots)
      {
        int & arm = pot.first;
        int & sta = pot.second;

        if(arm == 0)
          for(int pla = 2*nPlanes-1; pla >= 0; pla--)
          {
            int uv = pla % 2;

            bool hasDet = (hitDet[arm][sta][rpt][pla][reg] != 0);
            bool hasTra = (hitTra[arm][sta][rpt][pla][reg] != 0);

            print(hasDet,hasTra, uv);
          }
        else
          for(int pla = 0; pla < 2*nPlanes; pla++)
          {
            int uv = pla % 2;

            bool hasDet = (hitDet[arm][sta][rpt][pla][reg] != 0);
            bool hasTra = (hitTra[arm][sta][rpt][pla][reg] != 0);
            print(hasDet,hasTra, uv);
          }

        if(pot == pair<int,int>(0,0)) cerr << " |    | ";
                                 else cerr << " | ";
      }

      cerr << endl;
    }

    if(rpt == 0)
      cerr << " " << dots << "    " << dots << endl;
  }
  cerr << " " << dash << "    " << dash << endl;
}



