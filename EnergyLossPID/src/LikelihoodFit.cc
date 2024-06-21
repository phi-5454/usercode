#include "../interface/LikelihoodFit.h"
#include "../interface/Measurement.h"

#include <iostream>
#include <cstdio>
#include <cmath>

#define Sqr(x) ((x) * (x))

using namespace std;
using namespace CLHEP;

/*****************************************************************************/
inline double LikelihoodFit::sigmaD(double y)
{
  return sigma0 + b * y;
}

/*****************************************************************************/
double LikelihoodFit::probability
  (const Measurement & m, double epsilon)
{
  double y  = m.energy; // = g * y
  double sn = m.energySigma;
  double ls = m.path;   
  
  double Delta = m.epsilon * ls;
  
  double sD  = sigmaD(y);
  double sDD = sigmaD(Delta);

  double s   = sqrt(sD *sD  + sn*sn);
  double s0  = sqrt(sDD*sDD + sn*sn);
  
  double chi2;
  if(Delta < y - nu * s) chi2 = - 2*nu * (Delta - y)/s - nu*nu;
                    else chi2 =      Sqr((Delta - y))/ (s*s);

  return exp(-chi2/2) / (s / s0);
}

/*****************************************************************************/
void LikelihoodFit::calculate_wrt_y
  (const Measurement & m, double * val)
{
  double y  = m.energy; // = g * y
  double sn = m.energySigma;
  double ls = m.path;

  double Delta = m.epsilon * ls;

  double sD  = sigmaD(y);
  double sDD = sigmaD(Delta);
  double s   = sqrt(sD*sD + sn*sn);

  double d1 = (sDD*sD + sn*sn)/(s*s*s);
  double d2 = b*(sDD * (-2*sD*sD + sn*sn) - 3*sD * sn*sn)/(s*s*s*s*s);

  // Calculate derivatives with respect to y
  if(! m.isSaturated)
  { // normal
    if(Delta < y - nu * s)
    { // normal
      val[0] = - 2*nu * (Delta - y)/s - nu*nu;
      val[1] = 2*nu*d1;
      val[2] = 2*nu*d2;
    }
    else
    {
      val[0] = Sqr((Delta - y))/ (s*s);
      val[1] = 2*(-(Delta-y)/s)*d1;
      val[2] = 2*d1*d1 + 2*(-(Delta-y)/s)*d2;
    }
  }
  else
  { // saturated
    if(Delta < y + s)
    {
      val[0] = - (Delta - y)/s + 1;
      val[1] = d1;
      val[2] = d2;
    }
    else
    {
      val[0] = 0;
      val[1] = 0;
      val[2] = 0;
    }
  }
}

/*****************************************************************************/
void LikelihoodFit::calculate_wrt_epsilon
  (const Measurement & m, double epsilon, double * val)
{
  double y  = m.energy; // = g * y
  double sn = m.energySigma;
  double ls = m.path;

  double Delta = epsilon * ls;

  double sD  = sigmaD(y);

  double s   = sqrt(sD*sD + sn*sn);

  // Calculate derivatives with respect to Delta
  if(! m.isSaturated)
  { // normal
    if(Delta < y - nu * s)
    {
      val[0] = - 2*nu * (Delta - y)/s - nu*nu;
      val[1] = - 2*nu/s;
      val[2] = 0;
    }
    else
    {
      val[0] = Sqr((Delta - y))/ (s*s);
      val[1] = 2 * (Delta - y) / (s*s);
      val[2] =               2 / (s*s);
    }
  }
  else
  { // saturated
    if(Delta < y + s)
    {
      val[0] = - (Delta - y)/s + 1;
      val[1] = - 1/s;
      val[2] = 0;
    }
    else
    {
      val[0] = 0;
      val[1] = 0;
      val[2] = 0;
    }
  }

  // d/d Delta -> d/d epsilon
  val[1] *= ls;
  val[2] *= ls*ls;
}

