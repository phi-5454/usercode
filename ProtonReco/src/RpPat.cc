#include "../interface/RpPat.h"

#include "../interface/Parameters.h"

#include "../interface/RpFit.h"
#include "../interface/Helper.h"
#include "../interface/gzstream.h"

#include <algorithm>

using namespace std;

/*****************************************************************************/
// flag = 0 : collect and write
// flag = 1 : read patterns, above 0 permille
// flag = 2 : read patterns, dominant (>= 100)
// flag = 3 : read patterns, except 1nTu3 hit
// flag = 4 : don't read
RpPat::RpPat(const int flag, const int nCmTra) :
                  flag(flag), nCmTra(nCmTra)
{
  sFound = 0; sEmpty = 0;

  if(flag == 1)
  {
    readPatterns(0); // find highest num if pattern has no hole
    readPatterns(1); // take only those with num > highest num * minFraction
  }

  if(flag == 2 || flag == 3)
    readPatterns(flag);
}

/*****************************************************************************/
RpPat::~RpPat()
{
  if(flag == 0)
  {
    writePatterns();
    writeOccupancy();
  }
}

/*****************************************************************************/
void RpPat::writePatterns()
{
  cerr << Helper::col(5) << " writing roman pot patterns.."
       << Helper::col();

  char name[256];
  sprintf(name,"../out/rp/%dpart/patterns.out.gz",nCmTra);
  ogzstream file(name);

  for(auto & m : patterns)
  { 
    const auto & det   = m.first;
    const auto & vpats = m.second;

    file << " " << det.print() << endl;

    for(auto & v : vpats)
    {
      const auto & vpat = v.first;
      const auto & num  = v.second;

      for(auto & diff : vpat) // write all patterns!
        file << " " << diff;

      file << " " << num << endl;
    }
  }

  file.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
void RpPat::writeOccupancy()
{
  cerr << Helper::col(5) << " writing roman pot occupancy.."
       << Helper::col();

  char name[256];
  sprintf(name,"../out/rp/%dpart/occupancy.out.gz",nCmTra);
  ogzstream file(name);

  for(auto & m : occupancy)
  {
    const auto & det = m.first;

    file << " " << det.print() << endl;

    for(double strip = 0; strip < nStrips; strip += 0.5)
    {
      file << " " << strip;
 
      for(int pla = 0; pla < nPlanes; pla++)
        file  << " " << occupancy[det][pla][strip];

      file << endl;
    }

    file << endl << endl;
  }

  file.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
void RpPat::getDet(const string & line, RpDet & det)
{
  const vector<string> arg = Helper::parse(line,"|");

  int i = 0;
  det.arm = stoi(arg[i++]);
  det.sta = stoi(arg[i++]);
  det.rpt = stoi(arg[i++]);

  if(arg[i] == "u0") det.uv = 0;
                else det.uv = 1;
}

/*****************************************************************************/
void RpPat::readPatterns(int flag)
{
  cerr << Helper::col(2) << " reading roman pot patterns (flag "<<flag<<", "
                                              <<nCmTra<<" tracks).."
       << Helper::col();

  char name[256];
  sprintf(name,"../pars/patterns.out.gz");
  igzstream file(name);

  int nloc = 0, npat = 0;
  string line;
  RpDet det;

  while(getline(file, line))
  {
    const vector<string> a = Helper::parse(line," ");

    if(a.size() == 1) 
    {
      getDet(line, det);

      nloc++;
    }
    else 
    {
      const vector<string> arg = Helper::parse(line," ");
      vector<double> vpat;

      bool hasHole = false;
      for(size_t i = 0; i < arg.size()-1; i++)
      {
        double diff = stof(arg[i]);
        if(diff == Empty) hasHole = true;
        vpat.push_back(diff);
      }

      int num = stoi(arg.back());

      if(flag == 0 && !hasHole && num > highest[det]) // find highest
      { highest[det] = num; }

      if( (flag == 1 && num > highest[det] * 0e-3) || // above 0 permille
          (flag == 2 && num >= 1e+2) ||               // dominant
          (flag == 3 && !(det.print() == "0|0|0|u0"   // no 1nTu3 hit
                           && isFilled(vpat[2]))) )
      { patterns[det][vpat] = num; npat++; }
    } 
  }

  file.close();

  cerr << " [done " << nloc << " pots, " << npat << " patterns]" << endl;
}

/*****************************************************************************/
bool RpPat::isFilled(double x)
{
  return (x != Empty);
}

/*****************************************************************************/
int RpPat::nFilled(const double strip[])
{
  int nFilled = 0;
  for(int pla = 0; pla < nPlanes; pla++)
    if(isFilled(strip[pla])) nFilled++;

  return nFilled;
}

/*****************************************************************************/
bool RpPat::isSingle(double x)
{
  return (x == int(x));
}

/*****************************************************************************/
vector<double> RpPat::rebase(int ibas, const double clus[])
{
  // get base strip
  int sas = int(floor(clus[ibas]));

  vector<double> vpat;

  for(int pla = 0; pla < nPlanes; pla++)
  {
    const double & str = clus[pla];
    
    if(isFilled(str))
    { 
      double dif = str - sas; // rebase tracklet
      vpat.push_back(dif);
    }
    else
      vpat.push_back(Empty);
  }

  return vpat;
}

//
vector<double> RpPat::rebase(int ibas, const vector<double> & clus)
{
  // get base strip
  int sas = int(floor(clus[ibas]));

  vector<double> vpat;

  for(int pla = 0; pla < nPlanes; pla++)
  {
    const double & str = clus[pla];

    if(isFilled(str))
    {
      double dif = str - sas; // rebase tracklet
      vpat.push_back(dif);
    }
    else
      vpat.push_back(Empty);
  }

  return vpat;
}

/*****************************************************************************/
int RpPat::removeShear(int ibas, vector<double> & str)
{
  double sum = 0;
  int num = 0;

  for(int pla = ibas+1; pla < nPlanes; pla++) // no need for pla = ibas
  if(isFilled(str[pla]))
  { sum += str[pla]; num += pla-ibas; } // collect shear (weak mode) wrt ibas

  // take out shear (weak mode)
  int step = int(fabs(sum) / num) * (sum > 0 ? 1 : -1);

  for(int pla = ibas+1; pla < nPlanes; pla++)
  if(isFilled(str[pla]))
    str[pla] -= (pla-ibas)*step;

  return step;
}

/*****************************************************************************/
const map<RpDet, map<vector<double>,int>> & RpPat::getPatterns(bool addHoles)
{
  if(addHoles)
  {
    // look at all dets
    for(auto & m : patterns)
    {
      const auto & det  = m.first;
      const auto & vpat = m.second;

      // look at all patterns
      for(auto & v : vpat)
      {
        const auto & pat = v.first;
      //const auto & num = v.second;

        // is there a hole already?
        int ihole = -1;
        for(int i = 0; i < nPlanes; i++)
        if(!isFilled(pat[i]))
          ihole = i;

        if(ihole != -1)
        {
          // has a hole at ihole already
          // put in one more hole 
          for(int i = 0; i < nPlanes; i++)
          if(i != ihole) // no hole yet there
          {
            vector<double> str = pat;
            str[i] = Empty;

            str = rebase(getBase(str),str);

            patterns[det][str] = 0; // code
          }
        }
        else
        {
          // put in one hole 
          for(int i = 0; i < nPlanes; i++)
          {
            vector<double> str = pat;
            str[i] = Empty; 

            str = rebase(getBase(str),str);

            patterns[det][str] = 0; // code
          } 
       
          // put in two holes 
          for(int i = 0;   i < nPlanes-1; i++)
          for(int j = i+1; j < nPlanes  ; j++)
          {
            vector<double> str = pat;
            str[i] = Empty; 
            str[j] = Empty; 

            str = rebase(getBase(str),str);

            patterns[det][str] = 0; // code
          }
        }
      }
    }
  }

  return patterns;
}

/*****************************************************************************/
void RpPat::collectPatterns(const vector<RpTrack> & rpTracks)
{
  for(auto & track : rpTracks)
  {
    RpDet det = track.det;

    for(int uv = 0; uv < 2; uv++)
    {
      det.uv = uv;

      for(int pla = 0; pla < nPlanes; pla++)
      {
        mtx_occup.lock();
        occupancy[det][pla][track.clus[uv][pla]]++;
        mtx_occup.unlock();
      }

      int nFill = nFilled(track.clus[uv]);
      bool isSpecial = (det.print() == "0|0|0|u0");

      if( nFill == nPlanes || // (all filled) or (isSpecial and pla=2 miss)
         (nFill == nPlanes-1 && isSpecial && !isFilled(track.clus[uv][2])) )
      {
        // rebase wrt pla=0
        vector<double> vpat = rebase(0,track.clus[uv]); // 0 is always present

        // take out shear (weak mode)
        removeShear(0,vpat);

        // fill patterns
        mtx_pats.lock();
        patterns[det][vpat]++;
        mtx_pats.unlock();
      }
    }
  }
}

/*****************************************************************************/
int RpPat::getBase(const vector<double> & strip) // return index
{
  // find first filled plane, base
  int ibas = -1;
  for(int pla = 0; pla < nPlanes; pla++)
  if(isFilled(strip[pla]))
  { ibas = pla; break; }

  return ibas;
}

/*****************************************************************************/
vector<int> RpPat::toPattern(vector<double> & strip)
{
  // get base
  int ibas = getBase(strip);
  int base = int(strip[ibas]);

  // rebase
  vector<double> rem = rebase(ibas,strip);

  // remove shear
  int step = removeShear(ibas,rem);

  strip = rem; // write back

  return {ibas,base,step};
}

/*****************************************************************************/
// get three-element subsets of planes
void RpPat::getSubsets(const vector<double> & str,
                        vector<vector<int>> & sub)
{
  // 3 out of 5
  const vector<vector<int>> sets =
   { {0,1,2}, {0,1,3}, {0,1,4}, {0,2,3}, {0,2,4},
     {0,3,4}, {1,2,3}, {1,2,4}, {1,3,4}, {2,3,4} };

  // collect if filled
  for(auto & s : sets)
  if(isFilled(str[s[0]]) &&
     isFilled(str[s[1]]) &&
     isFilled(str[s[2]]))
    sub.push_back(s);
}

/*****************************************************************************/
/*
vector<double> RpPat::reduce(const vector<double> & str)
{
  int ibas = getBase(str);
  vector<double> rem = rebase(ibas,str);
  removeShear(ibas,rem);

  return rem;
*/

int RpPat::reduce(vector<double> & str)
{
   // find first filled plane, base
  int ibas = -1, sas = -1;

  for(int pla = 0; pla < nPlanes; pla++)
  if(isFilled(str[pla]))
  { ibas = pla; sas = int(floor(str[pla])); break; }

  for(int pla = ibas; pla < nPlanes; pla++)
  if(isFilled(str[pla]))
    str[pla] -= sas; // rebase tracklet

  double sum = 0;
  int num = 0;

  for(int pla = ibas+1; pla < nPlanes; pla++) // no need for pla = ibas
  if(isFilled(str[pla]))
  { sum += str[pla]; num += pla-ibas; } // collect shear (weak mode) wrt ibas

  // take out shear (weak mode)
  int step = int(fabs(sum) / num) * (sum > 0 ? 1 : -1);

  for(int pla = ibas+1; pla < nPlanes; pla++)
  if(isFilled(str[pla]))
    str[pla] -= (pla-ibas)*step;

  return step;
}

/*****************************************************************************/
void RpPat::collectForEffic(RpDet det, const vector<double> & str,
                            RpFit * theRpFit,
                            map<RpDet,double> & nFound,
                            map<RpDet,double> & nEmpty)
{
  bool isSpecial = (det.print() == "0|0|0|u0");

  // 0 or 1 holes, or 2 holes if special
  int nHoles = 0;
  int theHole = -1;
  for(int pla = 0; pla < nPlanes; pla++)
  if(!isFilled(str[pla]))
  {
    nHoles++;

    if(!isSpecial || pla != 2) // find the relevant hole
      theHole = pla;
  }

  if(!(nHoles < 2 || (nHoles == 2 && isSpecial))) return;

  // get base, rebase, remove shear
  int ibas = getBase(str);
  vector<double> rem = rebase(ibas,str);
  int step = removeShear(ibas,rem);

  if(!theRpFit->hasPattern(det,rem)) return;

  const double ds = 2.;       // PARAMETER
  const double minEff = 0.90; // PARAMETER

  auto & pats = patterns[det];

  //
  if(nHoles == 0 || (nHoles == 1 && isSpecial))
  {
    // take out one plane, see if there is a dominant pattern (>minEff = 90%)
    for(int pla = 0; pla < nPlanes; pla++)
    if((nHoles == 0 && !isSpecial) ||
       (nHoles == 1 &&  isSpecial && pla != 2))
    {
      //
      double sum = 0;      // sum of counts
      map<double,int> num; // counts in patterns
      for(double s = str[pla]-ds; s <= str[pla]+ds; s += 0.5)
      { 
        vector<double> mod = str; mod[pla] = s; // copy and set
        reduce(mod);
        auto it = pats.find(mod);

        if(it != pats.end()) num[s] = it->second;
                        else num[s] = 0;

        sum += num[s];
      }

      //
      const double & s = str[pla];
      if(isSingle(s))
      if(num[s] > minEff * sum)
      {
        det.pla = pla; det.str = s;
        mtx_found.lock(); nFound[det]++; sFound++; mtx_found.unlock();
        det.pla = -1 ; det.str = -1; // back
      }
    }
  }

  //
  if((nHoles == 1 && !isSpecial) ||
     (nHoles == 2 &&  isSpecial))
  {
    const int & pla = theHole;

    vector<pair<double,double>> fit = theRpFit->getPredicted(det, rem);

    // get predicted
    double pre = round(fit[pla].first)
               + int(floor(str[ibas])) + (pla-ibas)*step;

    // sum 
    double sum = 0;
    map<double,int> num;
    for(double s = pre-ds; s <= pre+ds; s += 0.5)
    {
      vector<double> mod = str; mod[pla] = s; // copy and set
      reduce(mod);
      auto it = pats.find(mod);

      if(it != pats.end()) num[s] = it->second;
                      else num[s] = 0;

      sum += num[s];
    }

    for(double s = pre-ds; s <= pre+ds; s += 0.5)
    if(isSingle(s))
    if(num[s] > minEff * sum)
    {
      det.pla = pla; det.str = s;
      mtx_empty.lock(); nEmpty[det]++; sEmpty++; mtx_empty.unlock();
      det.pla = -1 ; det.str = -1; // back

      break;
    }
  }
}

