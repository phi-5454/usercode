#include "../interface/FitProbabilities.h"

#include <cmath>
#include <fstream>

#include "../interface/MostProbable.h"

#include "TF1.h"
#include "TH2F.h"
#include "TGraphErrors.h"

#include "Math/MinimizerOptions.h"
#include "TVirtualFitter.h"

#include "../interface/MostProbable.h"

using namespace std;

static TVirtualFitter * fitter;

const int    etaBins = 20;
const double etaMin  = -1.;
const double etaMax  =  1.;

const int     ptBins = 40;
const double  ptMin  = 0.;
const double  ptMax  = 2.;

const int K = 3;

pair<int,pair<int,int> > gkey;

// charge, ieta, ipt, (logde, siglogde)
map<pair<int,pair<int,int> >, vector<pair<double,double> > > points;

// charge, ieta, ipt, part, (amp,mean)
map<pair<pair<int,int>,pair<int,int> >, double> amp, mean;
// charge, ieta, ipt, (scale)
map<pair<pair<int,int>,int>, double> scale;

double polpar[4];

int pass,code;
extern int gErrorIgnoreLevel;

/*****************************************************************************/
void loglikelihood
  (int &npar, double *gin, double &f, double *par, int flag)
{
  vector<double> mu(K), P(K);

  for(int k = 0; k < K; k++)
    mu[k] = par[k];

  P[0] = cos(par[K]);
  P[1] = sin(par[K]) * cos(par[K+1]);
  P[2] = sin(par[K]) * sin(par[K+1]);

  for(int k = 0; k < K; k++)
   P[k] = P[k] * P[k];

  double rescale = par[K+2];

  double logP = 0.;

  for(unsigned int n = 0; n < points[gkey].size(); n++)
  {
    double p = 0.;

    double sig = rescale * points[gkey][n].second; // don't upscale

    for(int k = 0; k < K; k++)
    if(P[k] > 0)
    {
      double q = (points[gkey][n].first - mu[k]) / sig;

      p += P[k] * (exp(-0.5*q*q) + 5e-3); // FIXME
    }

    p /= sig;

    if(p > 0)
      logP += -2 *log(p);
  }

  f = logP;
}


/*****************************************************************************/
FitProbabilities::FitProbabilities
  (map<pair<int,pair<int,int> >, 
       vector<pair<double,double> > > & points_)
{
  points = points_;

  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
  gErrorIgnoreLevel = 1001;

   // Initialize minimizer
  fitter = TVirtualFitter::Fitter(0, 2*K);
  fitter->SetFCN(loglikelihood);

  mostProbable = new MostProbable();
}

/*****************************************************************************/
FitProbabilities::~FitProbabilities()
{
}

/*****************************************************************************/
void FitProbabilities::printPropeller(int i)
{
  char c[4] = {'|', '/','-','\\'};

  cerr << c[i % 4] << "\b";
}

/*****************************************************************************/
double FitProbabilities::runLikelihood(int K,
   const vector<bool> & isFix,
   TVectorD & mu, TVectorD & P, TVectorD & mue)
{
  double arglist[100];

  arglist[0] = -1;
  fitter->ExecuteCommand("SET PRIN" , arglist, 1);
  fitter->ExecuteCommand("SET NOW"  , arglist, 0);

  double pre = 0.01;
  double tol = 0.2;
  fitter->SetParameter(0, "pion", mu[0], pre, mu[0]-tol, mu[0]+tol);
  fitter->SetParameter(1, "kaon", mu[1], pre, mu[1]-tol, mu[1]+tol);
  fitter->SetParameter(2, "prot", mu[2], pre, mu[2]-tol, mu[2]+tol);

  for(int k = 0; k < K; k++)
    P[k] = sqrt(P[k]);

  double theta = atan2(sqrt(P[1]*P[1] + P[2]*P[2]), P[0]);
  double phi   = atan2(P[2],P[1]);

  fitter->SetParameter(3, "theta",  theta, 0.01, 0., M_PI/2);
  fitter->SetParameter(4, "phi"  ,  phi  , 0.01, 0., M_PI/2);

  fitter->SetParameter(5, "rescale", 1., 0.01, 1./2, 2.);

  for(int k = 0; k < K; k++)
  if(isFix[k])
    fitter->FixParameter(k);

  arglist[0] = 0;
  fitter->ExecuteCommand("MIGRAD", arglist, 0);

  {
    double amin,edm,errdef;
    int nvpar,nparx;

    code = fitter->GetStats(amin, edm, errdef, nvpar, nparx);
  }

  for(int k = 0; k < K; k++)
  {
    mu[k]  = fitter->GetParameter(k);
    mue[k] = fitter->GetParError(k);
  }

  P[0] = cos(fitter->GetParameter(K));
  P[1] = sin(fitter->GetParameter(K)) * cos(fitter->GetParameter(K+1));
  P[2] = sin(fitter->GetParameter(K)) * sin(fitter->GetParameter(K+1));

  for(int k = 0; k < K; k++)
   P[k] = P[k] * P[k];

  return fitter->GetParameter(K+2);
}

/*****************************************************************************/
double polynom3(double *bg, double *a) // x and par
{
  double x = log(*bg);

  return exp(a[0] + a[1]*x + a[2]*x*x + a[3]*x*x*x);
}

/*****************************************************************************/
double FitProbabilities::getValue(double bg) // x and par
{
  return polynom3(&bg, &polpar[0]);
}

/*****************************************************************************/
void FitProbabilities::runLogLikelihood
  (const char* fileName,
   vector<double> & x,
   vector<double> & y,
   vector<double> & ye)
{
  const double mass[K] = { 0.139570, 0.493677, 0.938272 }; // GeV

  char name[256];
  sprintf(name,"../out/%s_%d.dat", fileName, pass);
  ofstream file(name);

  cerr << " log-likelihood";

  for(int icha = 0; icha < 2; icha++)
  {
    for(int ieta = 0; ieta < etaBins; ieta++)
    {
      cerr << ".";

      for(int ipt  = 0; ipt  <  ptBins; ipt++)
      {
        printPropeller(ipt);

        double eta = etaMin + (ieta + 0.5)/etaBins * (etaMax - etaMin);
        double pt  =  ptMin + (ipt  + 0.5)/ ptBins * ( ptMax -  ptMin);
        double p = pt * cosh(eta);

        pair<int,pair<int,int> > key(icha,pair<int,int>(ieta,ipt));

        if(points[key].size() > 1e+3) // need at least
        {
          file << " " << icha << " " << ieta << " " << ipt << " " << p;

          TVectorD mu(K);
          TVectorD mue(K);

          // set starting values
          for(int i = 0; i < K ; i++)
            if(pass == 1 || p/mass[i] < 0.2/mass[0]
                         || p/mass[i] > 1.9/mass[0])
              mu[i] = log(mostProbable->value(p/mass[i]));
            else
              mu[i] = log(getValue(p/mass[i]));

          TVectorD P(K);
          P[0] = 1 - 0.20 * p;
          P[1] =     0.12 * p;
          P[2] =     0.08 * p;

          gkey = key;

          // pass = 1
          vector<bool> isFix(K, false);

          if(pass == 2)
          {
            if(p > 1.10) isFix[1] = true; // fix kaon
            if(p > 1.70) isFix[2] = true; // fix prot
          }

          double rescale = runLikelihood(K, isFix, mu,P, mue);

          for(int i = 0; i < K ; i++)
          {
            file << " " << mu[i] << " " << mue[i] << " " << P[i];

            if( (i == 0 && p < 1.0 && p/mass[i] > 1.5) ||
                (i == 1 && p < 1.0 && p/mass[i] > 1.0) ||
                (i == 2 && p < 1.7 && p/mass[i] > 1.0) )
            if(code == 3 && mu[i] > log(2.5) && mu[i] < log(6.5) )
            {
               x.push_back(p/mass[i]);
               y.push_back(exp(mu[i]));
              ye.push_back(mue[i] * exp(mu[i]));
            }
          }

          file << " " << rescale
               << " " << points[key].size()
               << " " << code << endl;
        }
      }

      file << endl;
    }
    cerr << "|";

    file << endl;
  }
  file.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
void FitProbabilities::fitMostProbablePol3(const vector<double> & x,
                                      const vector<double> & y,
                                      const vector<double> & ye)
     // via a0,a1,a2,a3
{
  int n = 4;

  TF1 * func = new TF1("func", polynom3, -10.,10., n);

  int size = x.size();
  vector<double> xe(size, 0.0);
  vector<double> yee(size, 0.01);

  TGraphErrors * graph = new TGraphErrors(x.size(), &x[0],  &y[0],
                                                   &xe[0], &yee[0]);

  for(int i = 0; i < n; i++)
    func->SetParameter(i, polpar[i]);

  func->SetParName(0, "a0");
  func->SetParName(1, "a1");
  func->SetParName(2, "a2");
  func->SetParName(3, "a3");

  graph->Fit(func,"NR","");

  for(int i = 0; i < n; i++)
    cerr << " " << func->GetParameter(i)
         << " " << func->GetParError(i) << endl;

  // write back
  for(int i = 0; i < n; i++)
    polpar[i] = func->GetParameter(i);

  delete func;
  delete graph;
}

/*****************************************************************************/
void FitProbabilities::run()
{
  // Set initial parameters for pol3 fit
  polpar[0] =  1.50;
  polpar[1] = -0.85;
  polpar[2] =  0.45;
  polpar[3] = -0.07;

  pass = 1;

  // pass = 1
  {
    // Log-likelihood fit
    vector<double> x,y,ye;
    runLogLikelihood("mostprob", x,y,ye);

    // Fit most probable
    fitMostProbablePol3(x,y,ye);

    for(int i = 0; i < 4; i++)
      cerr << " pass=" << pass
           << " polpar[" << i << "] = " << polpar[i] << endl;

    pass++;
  }

  // pass = 2
  {
    // Log-likelihood fit
    vector<double> x,y,ye;
    runLogLikelihood("mostprob", x,y,ye);

    // Fit most probable
    fitMostProbablePol3(x,y,ye);

    for(int i = 0; i < 4; i++)
      cerr << " pass=" << pass
           << " polpar[" << i << "] = " << polpar[i] << endl;
  }
}

