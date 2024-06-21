#include "../interface/Fitters.h"

#include "../interface/LikelihoodFit.h"
#include "../interface/Measurement.h"
#include "../interface/SlimMeasurement.h"

#include <cmath>

#include <cstdio>
#include <fstream>

#define Sqr(x) ((x) * (x))

using namespace std;
using namespace CLHEP;

/*****************************************************************************/
Fitters::Fitters(const double sigma0_,
                 const double b_,
                 const vector<SlimMeasurement> & slim_) :
                 sigma0(sigma0_), b(b_), slim(slim_)
{
}

/*****************************************************************************/
Fitters::Fitters(const double sigma0_,
                 const double b_,
                 const vector<Measurement> & meas_) :
                 sigma0(sigma0_), b(b_), meas(meas_)
{
}

/*****************************************************************************/
Fitters::~Fitters()
{
}

/*****************************************************************************/
void Fitters::shft2(double &a, double &b, const double c)
{ a=b; b=c; }

/*****************************************************************************/
void Fitters::shft3(double &a, double &b, double &c, const double d)
{ a=b; b=c; c=d; }

/*****************************************************************************/
Measurement Fitters::toMeas(const SlimMeasurement & slim)
{
  Measurement meas;

  meas.energy      = slim.energy      * 5e-5;
  meas.energySigma = slim.energySigma * 5e-5;

  meas.path        = slim.path    * 1e-5;
  meas.isSaturated = slim.isSaturated;

  meas.epsilon     = slim.epsilon * 1e-3;

  meas.det = slim.det;

//  meas.nChannels = slim.nChannels;

  return meas;
}

/*****************************************************************************/
void Fitters::functionGain(double g, vector<double> & val)
{
  LikelihoodFit likelihoodFit(sigma0,b);

  for(int k = 0; k < 3; k++) val[k] = 0.;

  for(vector<SlimMeasurement>::const_iterator s = slim.begin();
                                              s!= slim.end(); s++)
  {
    double v[3];

    Measurement m = toMeas(*s);
    m.energy *= g;

    // d/d y
    likelihoodFit.calculate_wrt_y(m, v);

    double y = m.energy;

    // d/d g
    v[1] *=     y/g;
    v[2] *= Sqr(y/g);

    for(int k = 0; k < 3; k++) val[k] += v[k];
  }
}

/*****************************************************************************/
void Fitters::functionEpsilon(double epsilon, vector<double> & val)
{
  LikelihoodFit likelihoodFit(sigma0,b);

  for(int k = 0; k < 3; k++) val[k] = 0.;

//cerr << " --------- meas " << meas.size() << endl;

  for(vector<Measurement>::const_iterator m = meas.begin();
                                          m!= meas.end(); m++)
  {
    double v[3];
    likelihoodFit.calculate_wrt_epsilon(*m, epsilon, v);

    for(int k = 0; k < 3; k++) val[k] += v[k];
  }
}

/*****************************************************************************/
double Fitters::functionGain(double gain)
{
  vector<double> val(3);

  functionGain(gain, val);

  return val[0];
}

/*****************************************************************************/
double Fitters::scanModule(/*double function(double)*/)
{
  double fmin = 1e+30, xmin = 0.;

  for(double x = 0.01; x <= 100.; x *=1.1)
  {
    double f = functionGain(x); 

    if(f < fmin)
    { fmin = f; xmin = x; }
  }

  return xmin;
}

/*****************************************************************************/
double Fitters::goldenSearch(/*double function(double),*/
                             double ax, double bx, double cx)
{
  double xmin;

  const double R = 0.61803399; // (sqrt(5.) - 1) / 2
  const double C = 1 - R;

  double x1,x2;
  double x0=ax;
  double x3=cx;

  if(fabs(cx-bx) > fabs(bx-ax))
  { x1 = bx; x2 = bx+C*(cx-bx); }
  else
  { x2 = bx; x1 = bx-C*(bx-ax); }

  double f1 = functionGain(x1);
  double f2 = functionGain(x2);

  while(fabs(x3-x0) > 1e-3)
  {
    if(f2 < f1)
    {
      shft3(x0,x1,x2, R*x2+C*x3);
      shft2(   f1,f2, functionGain(x2));
    }
    else
    {
      shft3(x3,x2,x1, R*x1+C*x0);
      shft2(   f2,f1, functionGain(x1));
    }
  }

  if(f1 < f2) xmin=x1;
         else xmin=x2;

  return xmin;
}

/*****************************************************************************/
double Fitters::newtonMethodGain(pair<double,double> & value)
{
  int nStep = 0;

  double par = value.first; // input MeV/cm
  double dpar;

  vector<double> val(3);

  do
  {
    functionGain(par, val);

    if(val[2] != 0.) dpar = - val[1]/fabs(val[2]);
                else dpar = 1.;                     // step up, for epsilon

    if(par + dpar > 0) par += dpar; // ok
                  else par /= 2;    // half

    nStep++;
  }
  while(fabs(dpar) > 1e-3 && nStep < 50);

  value.first  = par;
  value.second = 2/val[2]; // sigma2 was 2

  return val[0];
}

/*****************************************************************************/
double Fitters::minimizeAllSaturated(pair<double,double> & value, int & nStep)
{
  nStep = 0;

  double par = 3.0; // input MeV/cm

  vector<double> val(3);

  do
  {
    functionEpsilon(par, val);


    if(val[1] != 0)
      par += - val[0]/val[1];
 
    nStep++;
  }
  while(val[0] > 1e-3 && val[1] != 0 && nStep < 10);

  value.first  = par;
  value.second = 0.; // CHECK

  value.second  = Sqr(value.first * 0.1); // 
  value.first  *= 1.1;                    //

  return val[0];
}

/*****************************************************************************/
double Fitters::newtonMethodEpsilon(pair<double,double> & value, int & nStep)
{
  nStep = 0;

  double par = 3.0; // input MeV/cm
  double dpar;

  vector<double> val(3);

  do
  {
    functionEpsilon(par, val);

    if(val[2] != 0.) dpar = - val[1]/fabs(val[2]);
                else dpar = 1.;                     // step up, for epsilon
                                                   
    if(par + dpar > 0) par += dpar; // ok
                  else par /= 2;    // half

    nStep++;
  }
  while(fabs(dpar) > 1e-3 && nStep < 50);

  value.first  = par;
  value.second = 2/val[2]; // sigma2 was 2!!

  return val[0];
}

