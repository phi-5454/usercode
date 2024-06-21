#include "../interface/ParticleType.h"
#include "../interface/MostProbable.h"

#include <iostream>
using namespace std;

#include "../interface/Enumerators.h"
#include "../interface/Masses.h"

/*****************************************************************************/
ParticleType::ParticleType(MostProbable * mostProbable_, int pass_) :
                             mostProbable(mostProbable_), pass(pass_)
{
}

/*****************************************************************************/
ParticleType::~ParticleType()
{
}

/*****************************************************************************/
// for depositMap and calibGain
int ParticleType::guess(double p, double y)
{
  int pid = unknown;

  double eps[nParts];

  double ylim;
  if(pass < 1) ylim = 3.20;
          else ylim = 3.05;

  for(int i = elec; i <= prot; i++)
    eps[i] = mostProbable->value(p/mass[i]);

  if(    y < (eps[elec] + eps[pion])/2 && p < 0.16)    pid = elec;
  else
  {
    if(  y < (eps[pion] + eps[kaon])/2 || y < ylim || p > 2)
    {
      if(p > mass[pion]) pid = pion; // only if p/m > 1
    }
    else
    {
      if(y > 3.7) // FIXME
      {
        if(y < (eps[kaon] + eps[prot])/2) { if(p < 0.70) pid = kaon; }
                                     else { if(p < 1.40) pid = prot; }
      }
    }
  }

  return pid;
}

/*****************************************************************************/
// for V0 selection
int ParticleType::sure (double p, double y)
{
  int pid = unknown;

  double eps[nParts];

  for(int i = elec; i <= prot; i++)
    eps[i] = mostProbable->value(p/mass[i]);

  if(y < (eps[pion] + eps[kaon])/2)
  { 
    if(p > 0.16 && p < 0.70) pid = pion;
  }
  else
  {
    if(y < (eps[kaon] + eps[prot])/2) { if(p < 0.70) pid = kaon; }
                                 else { if(p < 1.40) pid = prot; }
  }

  return pid;
}

