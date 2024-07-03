#include "../interface/SimplexMinimizer.h"

#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

/****************************************************************************/
SimplexMinimizer::SimplexMinimizer(int ndim) : ndim(ndim)
{
  mpts = ndim+1;
}

/****************************************************************************/
SimplexMinimizer::~SimplexMinimizer()
{
}

/****************************************************************************/
void SimplexMinimizer::get_psum(const Simplex & p, vector<double> & psum)
{
  for(int j=0; j<ndim; j++)
  {
    double sum = 0.0;

    for(int i=0; i<mpts; i++)
      sum += p[i][j];

    psum[j] = sum;
  }
}

/****************************************************************************/
double SimplexMinimizer::amotry(Simplex & p,
                                vector<double> & y, vector<double> & psum,
                                const int ihi, const double fac,
                                function<double(vector<double>)> & func)
// Helper function: Extrapolates by a factor fac through the face of the
// simplex across from the high point, tries it, and replaces the high point if
// the new point is better.
{
  vector<double> ptry(ndim);
  double fac1 = (1.0-fac)/ndim;
  double fac2 = fac1-fac;

  for(int j=0; j<ndim; j++)
    ptry[j] = psum[j]*fac1 - p[ihi][j]*fac2;

  double ytry = func(ptry); // Evaluate the function at the trial point.

  if(ytry < y[ihi])
  { // If it’s better than the highest, then replace the highest.
    y[ihi]=ytry;
 
    for(int j=0; j<ndim; j++)
    {
      psum[j] += ptry[j]  -p[ihi][j];
      p[ihi][j] = ptry[j];
    }
  }

  return ytry;
}

/****************************************************************************/
pair<vector<double>,double> SimplexMinimizer::minimize(
  Simplex p,
  function<double(vector<double>)> & func)
// Most general interface: initial simplex specified by the matrix
// p[0..ndim][0..ndim-1].  Its ndim+1 rows are ndim-dimensional vectors that
// are the vertices of the starting simplex.
{
  const int NMAX = 1000; // Maximum allowed number of function evaluations.
  const double TINY = 1.0e-10;

  int ihi,ilo,inhi;

  vector<double> psum(ndim),pmin(ndim,0),x(ndim);

  y.resize(mpts);

  for(int i=0; i<mpts; i++)
  {
    for(int j=0; j<ndim; j++)
      x[j]=p[i][j];

    y[i]=func(x);
  }

  int nfunc = 0;
  get_psum(p, psum);

  while(1)
  {
    ilo = 0;
    // First we must determine which point is the highest (worst),
    // next-highest, and lowest (best), by looping over the points in the
    // simplex.
    ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);

    for(int i=0; i<mpts; i++)
    {
      if(y[i] <= y[ilo]) ilo=i;
      if(y[i] > y[ihi]) { inhi=ihi; ihi=i; }
      else if(y[i] > y[inhi] && i != ihi) inhi=i;
    }

    double rtol = 2.0*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo])+TINY);

    // Compute the fractional range from highest to lowest and return if
    // satisfactory.
    if(rtol < ftol)
    { 
/*
      // If returning, put best point and value in slot 0.
      swap(y[0],y[ilo]);

      for(int i=0; i<ndim; i++)
      {
        swap(p[0][i], p[ilo][i]);
        pmin[i] = p[0][i];
      }

      const double & fmin = y[0];
*/
      // If returning, give average
      double fmin = 0;

      for(int i=0; i<mpts; i++)
      {
        fmin += y[i] / mpts;

        for(int j=0; j<ndim; j++)
          pmin[j] += p[i][j] / mpts;
      }

      return pair<vector<double>,double>(pmin,fmin);
    }

    if(nfunc >= NMAX)
    { 
//    cerr << " SimplexMinimizer: too many calls" << endl;

      vector<double> pmin(2); pmin = {-99,-99};
      double fmin = -99;

      return pair<vector<double>,double>(pmin,fmin);
    }

    nfunc += 2;

    // Begin a new iteration. First extrapolate by a factor 1 through the face
    // of the simplex across from the high point, i.e., reflect the simplex
    // from the high point.

    double ytry = amotry(p,y,psum,ihi,-1.0,func);

    if(ytry <= y[ilo])
    {
      // Gives a result better than the best point, so try an additional
      // extrapolation by a factor 2.
      ytry = amotry(p,y,psum,ihi,2.0,func);
    }
    else if(ytry >= y[inhi])
    {
      // The reflected point is worse than the second-highest, so look for an
      // intermediate lower point, i.e., do a one-dimensional contraction.
      double ysave = y[ihi];
      ytry = amotry(p,y,psum,ihi,0.5,func);

      if(ytry >= ysave)
      {
        // Can’t seem to get rid of that high point.
        for(int i=0; i<mpts; i++)
        {
          // Better contract around the lowest (best) point.
          if(i != ilo)
          {
            for(int j=0; j<ndim; j++)
              p[i][j] = psum[j] = 0.5*(p[i][j]+p[ilo][j]);

            y[i] = func(psum);
          }
        }

        nfunc += ndim;    // Keep track of function evaluations.
        get_psum(p,psum); // Recompute psum.
      }
    }
    else --nfunc; // Correct the evaluation count.
  }
  // Go back for the test of doneness and the next iteration.
}

/****************************************************************************/
pair<vector<double>,double> SimplexMinimizer::minimize
  (vector<double> point,
   vector<double> dels,
   function<double(vector<double>)> & func)
// Alternative interface that takes different displacements dels[0..ndim-1] in
// different directions for the initial simplex.
{
  Simplex p;

  vector<vector<double>> a;

  if(ndim == 2)
    a = { {-1./2, 1}, {-1./2,-1}, {1,0} };
  
  if(ndim == 3)
    a = { {1,1,1}, {1,-1,-1}, {-1,1,-1}, {-1,-1,1} };

  for(int i=0; i<ndim+1; i++)
  {
    p.push_back(point);

    for(int j = 0; j < ndim; j++)
      p[i][j] += a[i][j] * dels[j];
  }

  return minimize(p,func);
}

