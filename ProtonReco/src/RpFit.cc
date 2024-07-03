#include "../interface/RpFit.h"

#include "../interface/gzstream.h"
#include "../interface/Helper.h"
#include "../interface/SimplexMinimizer.h"
#include "../interface/RpPolygon.h"

#include "../interface/Parameters.h"

#include <cmath>
#include <iostream>
#include <iomanip>

#define sqr(x) ((x)*(x))

using namespace std;

/*****************************************************************************/
RpFit::RpFit(const map<RpDet, map<vector<double>,int>> & patterns,
             const int nCmTra, int flag) :
             patterns(patterns), nCmTra(nCmTra)
{
  //
  if(flag == 1)
    readShifts();

  if(flag == 2)
  {
    readShifts();
    readPatternFits();
  }
}

/*****************************************************************************/
RpFit::~RpFit()
{
}

/*****************************************************************************/
void RpFit::getDet(const string & line, RpDet & det)
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
void RpFit::readShifts()
{
  cerr << Helper::col(2)
       << " reading roman pot plane shifts (" << nCmTra << " tracks).."
       << Helper::col();

  string ntra = to_string(nCmTra);
  ifstream file("../pars/planeShifts_chi2.dat");

  int nDet = 0;
  string line;

  while(getline(file,line))
  {
    const vector<string> a = Helper::parse(line," ");

    RpDet det;
    getDet(a[0], det);

    for(size_t i = 1; i < a.size(); i++)
      shifts[det].push_back(stof(a[i]));

    nDet++;
  }

  file.close();

  cerr << " [done " << nDet << " pots]" << endl;
}

/*****************************************************************************/
bool RpFit::isFilled(double y) { return (y != Empty); }
bool RpFit::isSingle(double y) { return (y == int(y)); }

/*****************************************************************************/
vector<double> RpFit::getZ(const RpDet & det)
{
  // inside out
  vector<double> z;

  if(det.arm == 1)
  { // right arm
    if(det.uv == 0) z = {-20.3, -11.3, -2.3,  6.7, 15.7};
               else z = {-15.7,  -6.7,  2.3, 11.3, 20.3};
  }
  else
  { // left arm
    if(det.uv == 0) z = { 20.3,  11.3,  2.3, -6.7,-15.7};
               else z = { 15.7,   6.7, -2.3,-11.3,-20.3};
  }

  return z;
}

//
double RpFit::getDelta(double meas)
{
  const double fraction = 0.105; // from rp_double.gnu

  return (isSingle(meas) ? (1-fraction)/2 :
                              fraction /2); // in strip width units [pitch]
}

/*****************************************************************************/
double RpFit::getStrip(const RpDet & det, int pla, double a, double b)
{
  vector<double> shift = {0,shifts[det][0],shifts[det][1],shifts[det][2],0};

  vector<double> z = getZ(det);

  double pred = (a*z[pla] + b) + shift[pla];
  int istrip = int(pred + 0.5);

  const double fraction = 0.105; // from rp_double.gnu

  if(fabs(pred - istrip) < (1-fraction)/2)
    pred = istrip; // single
  else
  { // double
    if(pred > istrip) pred = istrip + 0.5;
                 else pred = istrip - 0.5;
  }

  return pred;
}

/*****************************************************************************/
double RpFit::localChi2(const vector<double> & pars)
{
  // refs
  const vector<double> & strip = my_strip;
  const vector<double> & shift = my_shift;

  // refs to pars
  const double & a = pars[0];
  const double & b = pars[1];

  //
  vector<double> z = getZ(my_det);

  // zero
  double val = 0;

  // collect chi2
  for(int pla = 0; pla < nPlanes; pla++)
  {
    const double & meas = strip[pla];

    if(isFilled(meas))
    {
      const double & x = z[pla];

      const double pred = (a*x + b) + shift[pla];

      double delta = getDelta(meas);

      double diff = 0.;

      if(pred < meas-delta) diff = pred - (meas-delta);
      else
      if(pred > meas+delta) diff = pred - (meas+delta);

      if(diff != 0.)
      {
        if(doSquare) val +=  sqr(diff);
                else val += fabs(diff);
      }
    }
  }

  return val;
}

/*****************************************************************************/
double RpFit::getTilt(const RpDet & det)
{
  double a;

  if(det.sta == 1) a = 45;   // far: 0|1 1|1 (was 220 far)
  else
  {
    if(det.arm == 0) a = 53; // near: 0|0    (was 210 far, tilted by 8 deg)
                else a = 37; // near: 1|0
  }

  return a * M_PI/180;
}

/*****************************************************************************/
Vector1 RpFit::stripToPos(const Vector1 & strip)
{
  const double no_of_strips_ = 512;
  const double pitch_ = 66E-3;
  const double y_width_ = 36.07;
  const double last_strip_to_border_dist_ = 1.4175;
  const double rec_u_0 = -14.9932;

  // GetHitPositionInReadoutDirection
  double position = last_strip_to_border_dist_ + (no_of_strips_-1)*pitch_
                    - y_width_/2. - strip.val * pitch_;

  // copy
  Vector1 pos;

  pos.val = position - rec_u_0;
  pos.sig = strip.sig * pitch_;

  return pos;
}

/*****************************************************************************/
Vector2 RpFit::getGlobalPosition(const RpDet & det,
                                 const vector<Vector1> & lpos)
{
  Vector1 u = stripToPos(lpos[0]); // u
  Vector1 v = stripToPos(lpos[1]); // v

  double a = getTilt(det);

  Vector2 gpos;

  int sign_x = (det.arm == 0 ? 1 : -1) * (det.sta == det.rpt ? 1 : -1);
  int sign_y = (det.rpt == 0 ? 1 : -1);

  //
  if(det.sta == 1)
  {
    gpos.x = (v.val*cos(a) - u.val*sin(a))*sign_x; // x
    gpos.y = (v.val*sin(a) + u.val*cos(a))*sign_y; // y

    gpos.cxx = sqr(v.sig*cos(a)) + sqr(u.sig*sin(a));
    gpos.cyy = sqr(v.sig*sin(a)) + sqr(u.sig*cos(a));

    gpos.cxy = (sqr(v.sig) - sqr(u.sig)) * cos(a)*sin(a) * sign_x*sign_y;
  }
  else
  {
    gpos.x = (u.val*cos(a) - v.val*sin(a))*sign_x; // x
    gpos.y = (u.val*sin(a) + v.val*cos(a))*sign_y; // y

    gpos.cxx = sqr(u.sig*cos(a)) + sqr(v.sig*sin(a));
    gpos.cyy = sqr(u.sig*sin(a)) + sqr(v.sig*cos(a));

    gpos.cxy = (sqr(u.sig) - sqr(v.sig)) * cos(a)*sin(a) * sign_x*sign_y;
  }

  return gpos;
}

/*****************************************************************************/
LocalFit RpFit::getLocalPosition(
  const RpDet & det, const vector<double> & strip, int ibas, int move, int step)
{
  if(fits[det].count(strip) > 0)
  {
    const vector<double> z = getZ(det);

    LocalFit localFit = fits[det][strip];

    // add back base and shear
    localFit.Cy += move + step*(0 - z[ibas])/(z[ibas+1] - z[ibas]);

    return localFit;
  }
  else
  {
#ifdef Debug
    cerr << " does not exist B :";

    for(auto & s : strip)
      cerr << " " << s;
#endif

    exit(1); 
  }
}

/*****************************************************************************/
vector<pair<double,double>> RpFit::getPredicted(
  const RpDet & det, const vector<double> & strip)//, int move, int step)
{
  vector<double> shift = {0,shifts[det][0],shifts[det][1],shifts[det][2],0};

  const vector<double> z = getZ(det);

  const LocalFit & f = fits[det][strip];

  const double & Cx = f.Cx;
  const double & Cy = f.Cy;
  const double & Vx = f.Vx;
  const double & Vy = f.Vy;

  vector<pair<double,double>> restored;

  for(int i = 0; i < nPlanes; i++) // unbase and add back shear
  {
    double mean  = Cy + Cx*z[i] + shift[i];  // mean
    double sigma = sqrt(Vy + Vx*sqr(z[i])); // sigma from Vy and Vx

    pair<double,double> v(mean,sigma);

    restored.push_back(v);
  }

  return restored;
}

/*****************************************************************************/
LocalFit RpFit::fitThroughPolygon(const vector<double> & strip,
                                        vector<vector<double>> & bands)
{
  const vector<double> z = getZ(my_det);
  const vector<double> & shift = my_shift;

  // set up bands
  bands.clear();

  for(int pla = 0; pla < nPlanes; pla++)
  if(isFilled(strip[pla]))
  {
    const double & meas = strip[pla];
    double delta = getDelta(meas);

    bands.push_back({-z[pla], meas - shift[pla], delta, double(pla)});
  }

  if(debug)
  {
    ofstream file("../out/rp/polygon/bands.dat");

    int i = 0;

    for(auto & band : bands)
    {
      const double & a = band[0];
      const double & b = band[1];
      const double & d = band[2];

      file << " " << a << " " << b+d << " "<< i << endl;
      file << " " << a << " " << b-d << " "<< i << endl;
      i++;
    }
    file.close();
  }

  //
  RpPolygon theRpPolygon;

  // Cx, Cy, sigma_x^2, sigma_y^2, value
  return theRpPolygon.findPolygon(bands, debug);
}

/*****************************************************************************/
double RpFit::getMean(const vector<double> & strip)
{
  // calculate mean
  double mean = 0;
  int n = 0;

  for(auto & s : strip)
  if(isFilled(s))
  { mean += s; n++; }

  mean /= n;

  return mean;
}

/*****************************************************************************/
vector<double> RpFit::fitThroughSimplex(const vector<double> & strip)
{
  SimplexMinimizer simplex(2);

  function<double(vector<double>)> funcLocalChi2 =
    bind(&RpFit::localChi2, this, std::placeholders::_1);

  // initial
  vector<double> point = {0,getMean(strip)};
  vector<double> dels  = {1,1};

  pair<vector<double>,double> result =
    simplex.minimize(point,dels, funcLocalChi2);

  vector<double> CxCyVal = {result.first[0], result.first[1], result.second};

  return CxCyVal;
}

/*****************************************************************************/
LocalFit RpFit::fitTracklet(const RpDet & det, vector<double> strip)
{
  my_strip = strip; 
  my_shift = {0,shifts[det][0],shifts[det][1],shifts[det][2],0};
  my_det   = det;

  // polygon
  vector<vector<double>> bands;
  LocalFit resPoly = fitThroughPolygon(strip, bands); // common area, get bands

  bool polyFail = (resPoly.Cx == 0. && resPoly.Cy == 0. &&
                   resPoly.Vx == 0. && resPoly.Vy == 0.);

  // simplex fit with sum ||
  doSquare = false;
  vector<double> resSimp = fitThroughSimplex(strip);
  doSquare = true;

  // no polygon, no common area
  if(polyFail)
  {
    // overwrite with simplex result
    resPoly.Cx  = resSimp[0]; // x
    resPoly.Cy  = resSimp[1]; // y
    resPoly.Vx  = 0;          //
    resPoly.Vy = sqr(0.30);   // PARAMETER
    resPoly.Val = resSimp[2]; // value
  }

  if(debug)
  {
    //
    ofstream file;

    //
    file.open("../out/rp/polygon/map.dat");

    vector<double> pars = {resSimp[0], resSimp[1]};
    vector<double> p(2);

    for(p[0] = pars[0]-0.1; p[0] < pars[0]+0.1; p[0] += 1e-3)
    {
    for(p[1] = pars[1]-1.; p[1] < pars[1]+1.; p[1] += 1e-2)
      file << " " << p[0] << " " << p[1] << " " << localChi2(p) << endl;
    file << endl;
    }

    file.close();

    //
    file.open("../out/rp/polygon/res.dat");

    file << " " << resPoly.Cx << " " << resPoly.Cy << " " << resPoly.Vy << endl;
    file << endl << endl;
    file << " " << resSimp[0] << " " << resSimp[1] << " " << resSimp[2] << endl;

    file.close();

    //
    // getchar();
    string command =
      "cd ../gnu ; ./gnuplot -e 'j="+to_string(iplot++)+"' rp_polygon.gnu";
    int result = system(command.c_str()); result++;
  }

  return resPoly;
}

/*****************************************************************************/
double RpFit::globalChi2(const vector<double> & shift)
{
  SimplexMinimizer simplex(2);

  function<double(vector<double>)> funcLocalChi2 =
    bind(&RpFit::localChi2, this, std::placeholders::_1);

  double gchi2 = 0;

  for(auto & m : patterns.at(my_det))
  {
    const auto & strip = m.first;
    const auto & num   = m.second;

    my_strip = strip; 
    my_shift = {0,shift[0],shift[1],shift[2],0};

    // initial
    vector<double> point = {0,getMean(strip)};
    vector<double> dels  = {1,1};

    pair<vector<double>,double> result =
      simplex.minimize(point,dels, funcLocalChi2);

//  const vector<double> & pars = result.first;
    const double         & val  = result.second;

    if(collectChi2)
      gchi2 += num * val; // collect local chi2
    else 
    if(val < 1e-6)
      gchi2 += -num;      // count number of zeros
  }

  return gchi2;
}

/*****************************************************************************/
void RpFit::optimizeShifts(bool collectChi2_)
{
  collectChi2 = collectChi2_;
  string turn = (collectChi2 ? "chi2" : "zero");
  string ntra = to_string(nCmTra);

  cerr << Helper::col(1)
       << " optimizing shifts ("
       << (collectChi2 ? "joint χ²" : "count zeros") << ")"
       << Helper::col() << endl;

  //
  SimplexMinimizer simplex(3);

  function<double(vector<double>)> funcGlobalChi2 =
    bind(&RpFit::globalChi2, this, std::placeholders::_1);

  // take all dets
  RpDet & det = my_det;

  ofstream     file("../out/rp/"+ntra+"part/planeShifts_"+turn+".dat");
  file << fixed << setprecision(3);
  cerr << fixed << setprecision(3);

  ofstream fileLine("../out/rp/"+ntra+"part/planeShifts_"+turn+"_line.dat");

  for(det.arm = 0; det.arm < 2; det.arm++)
  for(det.sta = 0; det.sta < 2; det.sta++)
  for(det.rpt = 0; det.rpt < 2; det.rpt++)
  for(det.uv  = 0; det.uv  < 2; det.uv++ )
  {
    file << " "  << det.print();
    cerr << "  " << det.print() << " :";

    vector<double> point = {0,0,0};
    vector<double> dels  = {1,1,1};

    pair<vector<double>,double> result =
      simplex.minimize(point,dels, funcGlobalChi2);

    const vector<double> & shift = result.first;
    const double         & gchi2 = result.second;

    // write and print
    for(int pla = 0; pla < nPlanes-2; pla++)
    {
      double rounded = int(shift[pla]*1000 + 0.5)/1000.;

      file << " " << setw(6) << rounded;
      cerr << " " << setw(6) << rounded;
    }
    file << endl;

    cerr << " (line..)";

    // line search, for demonstration
    for(int pla = 0; pla < nPlanes-2; pla++)
    {
      vector<double> pars = shift;

      for(pars[pla] = shift[pla] - 1;
          pars[pla] < shift[pla] + 1; pars[pla] += 5e-2)
      {
        double val = globalChi2(pars);
 
        if(fabs(val - gchi2) < 0.2 * fabs(gchi2))
          fileLine << " " << pars[pla] << " " << val << " " << pla+1 << endl;
      }

      fileLine << endl;
    }
    fileLine << endl;

    cerr << endl;
  }

  file.close();
  fileLine.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
bool RpFit::hasPattern(const RpDet & det, const vector<double> & str)
{
  return (fits[det].count(str) > 0);
}

/*****************************************************************************/
void RpFit::fitPatterns()
{
  cerr << Helper::col(1)
       << " fitting roman pot hit patterns.."
       << Helper::col() << endl;

  char name[256];
  sprintf(name,"../out/rp/%dpart/patterns.fit.gz",nCmTra);
  ogzstream file(name);

  srand48(1);

  for(auto & m : patterns)
  {
    const auto & det   = m.first;
    const auto & vpats = m.second;

    cerr << "  " << det.print();

    int n = 0;
    for(auto & v : vpats)
    {
      const auto & vpat = v.first;
      const auto & num  = v.second;

      debug = (nCmTra == 2 && drand48() < 1 && num > 8e+4 && iplot < 100);

      for(auto & v : vpat) // don't want emptry layer
        if(v == Empty) debug = false;

      LocalFit localFit = fitTracklet(det, vpat);
      debug = false;

      file << " " << det.print()
           << " " << num << " "
           << " " << localFit.Cx  // Cx
           << " " << localFit.Cy  // Cy
           << " " << localFit.Vx  // Vx
           << " " << localFit.Vy  // Vy  #6
           << " " << localFit.Val // Val #7
           << " ";

      for(auto & s : vpat)
        file << " " << s;

      file << endl;

      n++;
    }

    cerr << " [done " << n << " patterns]" << endl;
  }

  file.close();
}

/*****************************************************************************/
void RpFit::readPatternFits()
{
  cerr << Helper::col(2)
       << " reading roman pot fit patterns ("<<nCmTra<<" tracks).."
       << Helper::col();

  char name[256];
  sprintf(name,"../pars/patterns.fit.gz");
  igzstream file(name);

  int nFits = 0;
  string line;
  while(getline(file,line))
  {
    const vector<string> a = Helper::parse(line," ");

    RpDet det;
    getDet(a[0], det);
 
    const double & Cx  = stof(a[2]);
    const double & Cy  = stof(a[3]);
    const double & Vx  = stof(a[4]);
    const double & Vy  = stof(a[5]);
    const double & Val = stof(a[6]);

    vector<double> strip;
    for(int pla = 0; pla < nPlanes; pla++)
      strip.push_back(stof(a[7+pla]));

    fits[det][strip] = {Cx,Cy,Vx,Vy,Val};

    nFits++;
  }

  file.close();

  cerr << " [read " << nFits << "]" << endl;
}

