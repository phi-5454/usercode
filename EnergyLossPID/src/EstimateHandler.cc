#include "../interface/EstimateHandler.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <string>
#include <fstream>

#include "../interface/Measurement.h"
#include "../interface/SlimMeasurement.h"

#include "../interface/Fitters.h"
#include "../interface/Enumerators.h"

using namespace std;

/*****************************************************************************/
EstimateHandler::EstimateHandler()
{
  // set deposit parameters for all detector types
  sigma0 = 2e-3;
  b      = 0.095;
}

/*****************************************************************************/
int EstimateHandler::getDetId(uint32_t id, float thickness)
{
  int subdet = (id >> 25) & 0x7;

  subdet--;

  if(subdet == TEC3 && thickness > 400e-4)
    return TEC5;
  else
    return subdet;
}

/*****************************************************************************/
void EstimateHandler::calibrateGain
  (map<ChipId, vector<SlimMeasurement> > & hits,
   map<uint32_t, float> thickness, char * fileName,
   map<ChipId, float> & gain)
{
  ofstream fileGain(fileName);

  cerr << " Calibrate gain with golden section search + Newton method" << endl;
  cerr << " Number of keys : " << hits.size() << endl;
  cerr << "  ";

  int n = 0;

  for(map<ChipId, vector<SlimMeasurement> >::iterator
      key = hits.begin(); key!= hits.end(); key++)
  {
    if(n++ % (hits.size()/40) == 0) cerr << ".";

    if(key->second.size() >= 10) //
    if(key->first.first != 0)
    {
      uint32_t detId = (key->first).first;
      int det = getDetId(detId, thickness[detId]);
  
      // Initialize fitters, pass vector<SlimMeasurement>
      Fitters theFitters(sigma0, b, key->second);
  
      pair<double,double> value;
  
      // Golden section search [0.1,10.] with center 1.0
      value.first =  theFitters.goldenSearch(0.1, 1.0, 10.);
  
      // Newton method around value.first, get sigma2 = value.second
      theFitters.newtonMethodGain(value);

      float previousGain = (gain.count(key->first) > 0 ? gain[key->first] : 1);

      fileGain << " " << hex <<     key->first.first
               << " " << dec << int(key->first.second)
               << " " << dec <<     key->second.size()
               << " " << value.first  * previousGain
               << " " << value.second * previousGain*previousGain
               << " " << detNames[det] << endl;

      key->second.clear();
    }
  }

  cerr << endl << " [done]" << endl;
}

/*****************************************************************************/
bool EstimateHandler::allSaturated(const vector<Measurement> & meas)
{
  bool allSaturated = true;

  for(vector<Measurement>::const_iterator m = meas.begin();
                                          m!= meas.end(); m++)
    if(! m->isSaturated) allSaturated = false;

  return allSaturated;
}

/*****************************************************************************/
double EstimateHandler::estimate_(vector<Measurement> & meas,
                                  pair<double,double> & value)
{
  Fitters theFitters(sigma0,b, meas);

  int nStep;
  double chi2; 

  if(allSaturated(meas)) chi2 = theFitters.minimizeAllSaturated(value, nStep);
                    else chi2 = theFitters.newtonMethodEpsilon (value, nStep);

  return chi2;
}

/*****************************************************************************/
double EstimateHandler::estimate(vector<Measurement> & meas,
                                 pair<double,double> & value)
{
  int n = meas.size();

  // Compute original
  double chi2 = estimate_(meas,value);

  bool tryToRemove = ( n >= 3 &&
                       (chi2 > 1.3*n + 4*sqrt(1.3*n)) );

  // Try to remove lowest dE/dx measurement
  if(tryToRemove)
  {
    // Copy
    vector<Measurement> meas_(meas);
    vector<Measurement>::iterator drop;
  
    // Take lowest dE/dx hit
    double mindEdx = 1e+99;
    for(vector<Measurement>::iterator m = meas_.begin();
                                      m!= meas_.end(); m++)
    {
      double dEdx = m->energy / m->path;
  
      if(dEdx < mindEdx)
      { mindEdx = dEdx ; drop = m; } 
    }
  
    // Remove lowest dE/dx hit -> It will be absent during calibration!!
    meas_.erase(drop);
    pair<double,double> value_;
    double chi2_new  = estimate_(meas_,value_);
   
    // If got better, accept
    if(chi2_new < chi2 - 12) //
    {
      chi2 = chi2_new;
      meas = meas_;    // indeed remove the outlier!
      value = value_;
    } 
  }

  return chi2;
}

