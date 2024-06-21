#include "../interface/DataHandler.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "../interface/Measurement.h"
#include "../interface/SlimMeasurement.h"

#include "../interface/TPixelHit.h"
#include "../interface/TStripHit.h"

#include "../interface/TBunchCrossing.h"
#include "../interface/TVertex.h"
#include "../interface/TTrack.h"

#include "../interface/TPixelHit.h"
#include "../interface/TStripHit.h"

#include "../interface/MostProbable.h"
#include "../interface/ParticleType.h"
#include "../interface/EstimateHandler.h"

#include "../interface/FitStripCluster.h"

#include "../interface/FitProbabilities.h"

#include "TFile.h"
#include "TH2S.h"
#include "TH2F.h"
#include "TH3I.h"

using namespace std;

#include "../interface/Masses.h"

enum { isNormal, isBelow, isOver };

#define minHitsSub 2
#define minHitsAll 4

#define minHitsCal 6
#define minHitsFit 4

const int    lpBins = 150;
const double lpMin  = -3.;
const double lpMax  =  3.;

const int    ldeBins = 175;
const double ldeMin  =  0.0;
const double ldeMax  =  3.5;

const int    etaBins = 20;
const double etaMin  = -1.;
const double etaMax  =  1.;

const int     ptBins = 40;
const double  ptMin  = 0.;
const double  ptMax  = 2.;

//
const int    lbgBins =  100;
const double lbgMin  = -1.25;
const double lbgMax  =  2.75;

const int    pathBins = 150;
const double pathMin  =    0e-4;
const double pathMax  = 1500e-4;

const int    enerBins = 200;
const double enerMin  = 0.; // MeV
const double enerMax  = 1.;

const int    esigBins = 100;
const double esigMin  = 0.; // MeV
const double esigMax  = 0.1;


/*****************************************************************************/
DataHandler::DataHandler()
{
  pass = -99;

  coupledHits = false;
  calibGain   = false;
  depositMap  = false;
  applyGain   = false;
  surveyMap   = false;
  fitProbs    = false;

  mostProbable = new MostProbable();
  particleType = new ParticleType(mostProbable, pass);
}

/*****************************************************************************/
DataHandler::DataHandler(std::string tag_) : tag(tag_)
{
  pass = -99;

  coupledHits = false;
  calibGain   = false;
  depositMap  = false;
  applyGain   = false;
  surveyMap   = false;
  fitProbs    = false;

  mostProbable = new MostProbable();
  particleType = new ParticleType(mostProbable, pass);
}

/*****************************************************************************/
DataHandler::~DataHandler()
{
}

/*****************************************************************************/
void DataHandler::options(int arg, char **arc, int chunk)
{
  int i = 1;

  do
  {
    if(strcmp(arc[i],"-coupledHits") == 0) coupledHits = true;
    if(strcmp(arc[i],"-calibGain"  ) == 0) calibGain   = true;
    if(strcmp(arc[i],"-depositMap" ) == 0) depositMap  = true;
    if(strcmp(arc[i],"-applyGain"  ) == 0) applyGain   = true;
    if(strcmp(arc[i],"-surveyMap"  ) == 0) surveyMap   = true;
    if(strcmp(arc[i],"-fitProbs"   ) == 0) fitProbs    = true;

    if(strcmp(arc[i],"-pass" ) == 0) pass = atoi(arc[++i]);

    i++;
  }
  while(i < arg);

  if(chunk != 1) fitProbs = false;
}

/*****************************************************************************/
void DataHandler::beginJob()
{
  if(pass == -1)
    cerr << "\033[22;31m" << "Collecting coupled strip hits (pass -1)..."
         << "\033[22;0m"  << endl;

  if(pass == 0)
    cerr << "\033[22;31m" << "Collecting hits for gain calibration (pass 0)..."
         << "\033[22;0m"  << endl;

  if(pass == 1)
    cerr << "\033[22;31m" << "Collecting hits for gain calibration (pass 1)..."
         << "\033[22;0m"  << endl;

  if(pass == 2)
    cerr << "\033[22;31m" << "Estimating most probable energy loss (pass 2)..."
         << "\033[22;0m"  << endl;

  if(applyGain ) readGainCorrection();

  if(!coupledHits)
  {
    readStripProps();
    estimateHandler = new EstimateHandler();
  }

  if(coupledHits || depositMap || surveyMap)
  {
  cerr << " creating histos.. (";

  char name[256];

  //
  if(coupledHits)
  {
    cerr << "c";

    for(int det = 0; det < nDets; det++)
    {
      sprintf(name,"coupling_%s", detNames[det].c_str());

      TH2S h(name, name, 256, -0.5, 255.5,
                         256, -0.5, 255.5);
      hcoupling.push_back(h);
    }
  }

  if(depositMap)
  {
    cerr << "de";

    for(int cha  = pos; cha  < nCharges; cha++ )
    for(int part = 0  ; part < nParts  ; part++)
    for(int det  = 0  ; det  < nDets   ; det++ )
    {
      sprintf(name,"deposit_%s_%s_%s", chargeNames[cha].c_str(),
                                        partNames[part].c_str(),
                                          detNames[det].c_str());

      TH3I hd(name, name,  lbgBins, lbgMin, lbgMax,
                         pathBins, pathMin, pathMax,
                         enerBins, enerMin, enerMax);
      deposit3[cha].push_back(hd);

      sprintf(name,"esigma_%s_%s_%s", chargeNames[cha].c_str(),
                                       partNames[part].c_str(),
                                         detNames[det].c_str());

      TH3I he(name, name,  lbgBins, lbgMin, lbgMax,
                         pathBins, pathMin, pathMax,
                         esigBins, esigMin, esigMax);
      esigma3[cha].push_back(he);
    }
  }

  //
  if(surveyMap)
  {
    cerr << "s";

    //
    int j = 0;

    for(int ver  = 0; ver  < nVers ; ver++ )
    {
      for(int cha  = 0; cha  < nCharges; cha++ )
      {
        sprintf(name,"histos_%s_%s", verNames[ver].c_str(),
                                  chargeNames[cha].c_str());

        TH2F h(name,name, lpBins,  lpMin,  lpMax,
                         ldeBins, ldeMin, ldeMax);

        histos[j].push_back(h);
      }

      j++;
    }
  }

  cerr << ") [done]" << endl;
  }
}

/*****************************************************************************/
void DataHandler::writeHistos(const vector<TH2S> & histos)
{
  for(vector<TH2S>::const_iterator h2 = histos.begin();
                                   h2!= histos.end(); h2++)
    h2->Write();
}

/*****************************************************************************/
void DataHandler::writeHistos(const vector<TH2F> & histos)
{
  for(vector<TH2F>::const_iterator h2 = histos.begin();
                                   h2!= histos.end(); h2++)
  {
    h2->Write();
  }
}

/*****************************************************************************/
void DataHandler::writeHistos(const vector<TH3I> & histos)
{
  for(vector<TH3I>::const_iterator h3 = histos.begin();
                                   h3!= histos.end(); h3++)
  {
    h3->Write();
  }
}

/*****************************************************************************/
void DataHandler::endJob(int chunk)
{
  // Perform gain calibration
  if(calibGain)
  {
    char fileName[256];

    // Write hits to file
    sprintf(fileName,"../out/hits_%d.bin", chunk);

    ofstream outFile(fileName, ios::binary);

    // hits
    {
      int size = hits.size();
      outFile.write((char *)&size, sizeof(size));
  
      for(map<ChipId, vector<SlimMeasurement> >::const_iterator
          hit = hits.begin(); hit != hits.end(); ++hit)
      {
        outFile.write((char *)&hit->first , sizeof(hit->first ));

        int n = hit->second.size();
        outFile.write((char *)&n, sizeof(n));

        for(vector<SlimMeasurement>::const_iterator
          slim = hit->second.begin(); slim!= hit->second.end(); slim++)
        outFile.write((char *)&(*slim), sizeof(*slim));
      }
    }

    // thickness
    {
      int size = thickness.size();
      outFile.write((char *)&size, sizeof(size));

      for(map<uint32_t, float>::const_iterator
          thick = thickness.begin(); thick != thickness.end(); ++thick)
      {
        outFile.write((char *)&thick->first , sizeof(thick->first ));
        outFile.write((char *)&thick->second, sizeof(thick->second));
      }
    }

    outFile.close();
  }

  if(coupledHits || depositMap || surveyMap)
  {
  // Write out histos
  char fileName[256];
  sprintf(fileName,"../out/data_%d_%d.root",pass,chunk);

  TFile * file = new TFile(fileName,"recreate");

  cerr << " writing histos.. (";
 
  //
  if(coupledHits)
  {
    cerr << "c"; writeHistos(hcoupling);
  }

  // 
  if(depositMap)
  {
    cerr << "d";
    for(int i = 0; i < nCharges; i++) writeHistos(deposit3[i]);

    cerr << "e";
    for(int i = 0; i < nCharges; i++) writeHistos(esigma3[i]);
  }

  //
  if(surveyMap)
  {
    cerr << "s";
    for(int i = 0; i < nVers; i++) writeHistos(histos[i]);
  }

  cerr << ") [done]" << endl;

  //
  file->Close();
  }

  //
  if(fitProbs)
  {
    FitProbabilities fitProbabilities(points);

    fitProbabilities.run();
  }
}

/*****************************************************************************/
void DataHandler::readStripProps()
{
  cerr << " reading strip properties from";

  char fileName[256];

  if(pass != -99)
  {
    sprintf(fileName,"UserCode/EnergyLossPID/data/stripProps.par");
    cerr << " stripProps.par.." << endl;
  }
  else
  {
    sprintf(fileName,"UserCode/EnergyLossPID/data/stripProps_%s.par",
            tag.c_str());
    cerr << " stripProps_" << tag << ".par..";
  }

  edm::FileInPath fileInPath(fileName);
  ifstream file(fileInPath.fullPath().c_str());

  int det;
  for(det = PXB; det <= PXF; det++)
  {
      thr.push_back(0.);
    alpha.push_back(0.);
    sigma.push_back(0.);
  }
  
  while(!file.eof())
  {
    string detName;
    float f;

    file >> detName;
    file >> f;   thr.push_back(f);
    file >> f; alpha.push_back(f);
    file >> f; sigma.push_back(f);

    det++;
  }

  file.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
int DataHandler::getDetId(uint32_t id, float thickness)
{
  int subdet = (id >> 25) & 0x7;

  subdet--;

  if(subdet == TEC3 && thickness > 400e-4)
    return TEC5;
  else
    return subdet;
}

/*****************************************************************************/
SlimMeasurement DataHandler::toSlim(const Measurement & meas)
{
  SlimMeasurement slim;

  slim.energy      = int(meas.energy      * 2e+4);
  slim.energySigma = int(meas.energySigma * 2e+4);

  slim.path        = int(meas.path    * 1e+5);
  slim.isSaturated =     meas.isSaturated;

  slim.epsilon     = int(meas.epsilon * 1e+3);

  slim.det = meas.det;

  return slim;
}

/*****************************************************************************/
void DataHandler::readGainCorrection()
{
  cerr << " reading gain from";

  char fileName[256];

  if(pass != -99)
  {
    sprintf(fileName,"UserCode/EnergyLossPID/data/gain_%d.dat",pass-1);
    cerr << " gain_" << pass-1 << ".dat"; 
  }
  else
  {
    sprintf(fileName,"UserCode/EnergyLossPID/data/gain_%s.dat",tag.c_str());
    cerr << " gain_" << tag << ".dat"; 
  }

  edm::FileInPath fileInPath(fileName);
  ifstream fileGain(fileInPath.fullPath().c_str());

  int i = 0;
  while(!fileGain.eof())
  {
    uint32_t det;
    int chip;

    int d;
    float g,f;
    string s;

    fileGain >> hex >> det;
    fileGain >> dec >> chip;

    ChipId detId(det, (unsigned char)chip);

    fileGain >> dec >> d;
    fileGain >> g; fileGain >> f;
    fileGain >> s;

    if(!fileGain.eof())
    {
      if(g > 0.5 && g < 2.0) gain[detId] = g;
                        else gain[detId] = -1.;
    }

    if(i++ % 5000 == 0) cerr << ".";
  }

  fileGain.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
void DataHandler::correctEnergy(float & energy, const ChipId & detId)
{
  if(applyGain)
  if(gain.count(detId) > 0)
    energy *= gain[detId];
}

/*****************************************************************************/
double DataHandler::effectivePath(double l, int detType)
{
  double a  = 0.07;   //
  double l0 = 450e-4; // cm

  return l * (1 + a *log(l/l0));
}

/*****************************************************************************/
vector<Measurement> DataHandler::getPixelMeasurements(const TTrack & track)
{
  vector<Measurement> meas;

  for(vector<TPixelHit>::const_reverse_iterator
        hit = track.pixelHits.rbegin();
        hit!= track.pixelHits.rend(); hit++)
  {
    thickness[hit->detId] = hit->thickness;

    if(hit->forEloss)
    {
      Measurement m;

      m.energy = hit->Delta; 
      m.chipId = ChipId(hit->detId,hit->chip);

      if(pass < 2 || gain[m.chipId] > 0) 
      {
        correctEnergy(m.energy, m.chipId);

        m.det         = getDetId(hit->detId, hit->thickness);
        m.path        = effectivePath(hit->x, m.det);

        m.isSaturated = hit->isSaturated; // pixel, copy

        // Estimate sigma for cluster
        m.nChannels   = hit->nChannels;
        m.energySigma = 10e-3 * sqrt(hit->nChannels);

        m.r = hit->r;
        m.p = hit->p;

        if(m.det == PXF || m.energy > 50e-3) // THRESHOLD
          meas.push_back(m);
      }
    }
  }

  return meas;
}

/*****************************************************************************/
vector<Measurement> DataHandler::getStripMeasurements(const TTrack & track)
{
  vector<Measurement> meas;

  for(vector<TStripHit>::const_reverse_iterator
        hit = track.stripHits.rbegin();
        hit!= track.stripHits.rend(); hit++)
  {
    thickness[hit->detId] = hit->thickness;

    if(hit->forEloss && hit->meas <= int(fabs(hit->pred)) + 2 + 2)
    {
      Measurement m;
      m.det  = getDetId(hit->detId, hit->thickness);
      m.chipId = ChipId(hit->detId, hit->chip);

      if(pass < 2 || gain[m.chipId] > 0) 
      {
        m.path   = effectivePath(hit->x, m.det);

        // Set threshold
        double threshold = thr[m.det];

        // Fill measured strip deposits
        vector<pair<double,int> > b;

        b.push_back(pair<double,int>(threshold, isBelow));

        for(vector<uint8_t>::const_iterator adc = hit->adc.begin();
                                            adc!= hit->adc.end(); adc++)
        { 
          if(*adc > 253)
            b.push_back(pair<double,int>(254.              , isOver  ));
          else 
            b.push_back(pair<double,int>(double(*adc + 0.5), isNormal));
        }

        b.push_back(pair<double,int>(threshold, isBelow));

        // Fit
        FitStripCluster theFitter(b, alpha[m.det],sigma[m.det]);

        double chi2; pair<double,double> result;

        int width = theFitter.run(chi2, result);
        m.nChannels = width;

        // ADC -> e- -> MeV
        m.energy      =      result.first   * 217 * 3.61e-6;
        m.energySigma = sqrt(result.second) * 217 * 3.61e-6;

        m.isSaturated = false;  // strip, not saturated

        correctEnergy(m.energy     , m.chipId);
        correctEnergy(m.energySigma, m.chipId);

        m.r = hit->r;
        m.p = hit->p;

        meas.push_back(m);
      }
    }
  }

  return meas;
}

/*****************************************************************************/
void DataHandler::getCoupledMeasurements(const TTrack & track)
{
  for(vector<TStripHit>::const_iterator hit = track.stripHits.begin();
                                        hit!= track.stripHits.end(); hit++)
  if(hit->forEloss && hit->meas <= int(fabs(hit->pred)) + 2 + 2)
  if(fabs(hit->pred) < 0.1)
  {
    Measurement m;
    m.det  = getDetId(hit->detId, hit->thickness);

    if(hit->adc.size() == 2)
      hcoupling[m.det].Fill( max(hit->adc[0], hit->adc[1]),
                             min(hit->adc[0], hit->adc[1]) );

    if(hit->adc.size() == 3)
    {
      hcoupling[m.det].Fill(hit->adc[1], hit->adc[0]);
      hcoupling[m.det].Fill(hit->adc[1], hit->adc[2]);
    }
  }
}

/*****************************************************************************/
void DataHandler::fillHits(vector<Measurement> & meas,
                           int cha, double bg)
{
  for(vector<Measurement>::iterator m = meas.begin();
                                    m!= meas.end(); m++)
  {
    if(hits.count(m->chipId) == 0)
       hits[m->chipId].reserve(500);

    if(hits[m->chipId].size() < 1000) // FIXME
       hits[m->chipId].push_back(toSlim(*m));
  }
}

/*****************************************************************************/
void DataHandler::copyEpsilon(vector<Measurement> & meas,
                              int pid)
{
  int i = 0;
  for(vector<Measurement>::iterator m = meas.begin();
                                    m!= meas.end(); m++)
    m->epsilon = epsilon[pid][i++];
}

/*****************************************************************************/
void DataHandler::fillHistos(const TTrack & track,
                             const vector<int> & nhits,
                             const vector<pair<double,double> > & values,
                             const vector<Measurement> & measAll)
{
  double p = track.pt * cosh(track.eta);
  double logp = log(p);

  int cha = (track.charge > 0 ? pos : neg);

  // Deposit
  if(depositMap && nhits[all] >= minHitsAll)
  {
    int pid = particleType->guess(p,values[all].first);

    if(pid != unknown)
    {
      for(vector<Measurement>::const_iterator m = measAll.begin();
                                              m!= measAll.end(); m++)
      {
        double lbg = log(m->p / mass[pid]);

        deposit3[cha][pid*nDets + m->det].Fill(lbg, m->path, m->energy);
         esigma3[cha][pid*nDets + m->det].Fill(lbg, m->path, m->energySigma);

        deposit3[cha][unknown*nDets + m->det].Fill(lbg, m->path, m->energy);
         esigma3[cha][unknown*nDets + m->det].Fill(lbg, m->path, m->energySigma);
      }
    }
  }

  // Histos
  for(int ver = 0; ver < nVers; ver++) // pix str all
  if( (ver != all && nhits[ver] >= minHitsSub) ||
      (ver == all && nhits[ver] >= minHitsAll) )
  {
    double lde = log(values[ver].first);

    histos[ver][cha].Fill(logp, lde);
  }
}

/*****************************************************************************/
pair<double,double> DataHandler::processTrack(const TTrack & track,
                               vector<reco::DeDxData> & estimatePix,
                               vector<reco::DeDxData> & estimateStr,
                               vector<reco::DeDxData> & estimateAll)
{
  vector<Measurement> measPix, measStr, measAll;
  pair<double,double> valuePix(0.,0.), valueStr(0.,0.), valueAll(0.,0.);

  // Get hits, estimate
  measPix = getPixelMeasurements(track);
  measStr = getStripMeasurements(track);

  measAll.insert(measAll.end(), measPix.begin(), measPix.end());
  measAll.insert(measAll.end(), measStr.begin(), measStr.end());

  unsigned int nhits;

  // Pix
  nhits = measPix.size();
  if(measPix.size() >= minHitsSub)
    estimateHandler->estimate(measPix, valuePix);

  estimatePix.push_back(reco::DeDxData(valuePix.first,
                                       sqrt(valuePix.second), nhits));

  // Str
  nhits = measStr.size();
  if(measStr.size() >= minHitsSub)
    estimateHandler->estimate(measStr, valueStr);

  estimateStr.push_back(reco::DeDxData(valueStr.first,
                                       sqrt(valueStr.second), nhits));

  // All
  nhits = measAll.size();
  if(measAll.size() >= minHitsSub)
    estimateHandler->estimate(measAll, valueAll);

  estimateAll.push_back(reco::DeDxData(valueAll.first,
                                       sqrt(valueAll.second), nhits));

  return valueAll;
}

/*****************************************************************************/
void DataHandler::processTrack(const TTrack & track)
{
  if(coupledHits)
  {
    getCoupledMeasurements(track);
    return;
  }

  vector<Measurement> measPix, measStr, measAll;
  pair<double,double> valuePix(0.,0.), valueStr(0.,0.), valueAll(0.,0.);

  // Get hits, estimate
  measPix = getPixelMeasurements(track);
  measStr = getStripMeasurements(track);

  measAll.insert(measAll.end(), measPix.begin(), measPix.end());
  measAll.insert(measAll.end(), measStr.begin(), measStr.end());

  if(measPix.size() >= minHitsSub)
    estimateHandler->estimate(measPix, valuePix);
  if(measStr.size() >= minHitsSub)
    estimateHandler->estimate(measStr, valueStr);
  if(measAll.size() >= minHitsAll)
    estimateHandler->estimate(measAll, valueAll);

  if(depositMap || surveyMap || calibGain || fitProbs)
  {
    double p = track.pt * cosh(track.eta);

    // Fixed epsilon for particle types
    epsilon.clear(); epsilon.resize(nParts);
  
    for(unsigned int i = 0; i < measAll.size(); i++)
    {
      // Calculate hit epsilon with different mass assumptions
      for(int pid = elec; pid < nParts; pid++)
        epsilon[pid].push_back(mostProbable->value(measAll[i].p / mass[pid]));
    }
  
    //
    if(depositMap || surveyMap)
    {
      vector<int> nhits; vector<pair<double,double> > values;
  
      nhits.push_back(measPix.size()); values.push_back(valuePix);
      nhits.push_back(measStr.size()); values.push_back(valueStr);
      nhits.push_back(measAll.size()); values.push_back(valueAll);

      fillHistos(track, nhits, values, measAll);
    }
  
    //
    if(calibGain && measAll.size() >= minHitsCal)
    {
      int pid = particleType->guess(p,valueAll.first);
  
      if(pid != unknown)
      {
        copyEpsilon(measAll, pid);
  
        fillHits(measAll, (track.charge > 0 ? pos : neg), p/mass[pid]);
      }
    }

    //
    if(fitProbs)
    {
      double eps    = valueAll.first;
      double sigma2 = valueAll.second;
  
      // transform to logde and its sigma
      double    logde = log(eps);
      double siglogde = sqrt(sigma2)/eps;
  
      bool isOk = (!(std::isnan(   logde) || std::isinf(   logde) ||
                     std::isnan(siglogde) || std::isinf(siglogde)));
  
      if(measStr.size() > 0 && 
         measAll.size() >= minHitsFit &&
         isOk && logde > ldeMin && logde < ldeMax)
      {
        double eta = track.eta;
        double pt  = track.pt;

        int icha = (track.charge > 0 ? pos : neg);

        if(fabs(eta) < 1)
        {
          int ieta = int((eta - etaMin)/(etaMax - etaMin) * etaBins);
          int ipt  = int(( pt -  ptMin)/( ptMax -  ptMin) *  ptBins);
  
          pair<double,double> point(logde, siglogde);
          pair<int,pair<int,int> > key(icha, pair<int,int>(ieta,ipt));

          points[key].push_back(point);
        }
      }
    }
  }
}

/*****************************************************************************/
void DataHandler::processEvent(const TBunchCrossing * bunx)
{
  // Take at most bunx with two vertices
//  if(bunx->recVertices.size() <= 2) // FIXME
  for(vector<TVertex>::const_iterator vertex = bunx->recVertices.begin();
                                      vertex!= bunx->recVertices.end();
                                      vertex++)
  // Take only the main vertex
//  if(vertex == bunx->recVertices.begin()) FIXME
  {
    // Primary tracks
    for(vector<TTrack>::const_iterator
        track = vertex->tracks.begin();
        track!= vertex->tracks.end(); track++)
//    if(/*track->isPrimary &&*/ track->isHighPurity)
      processTrack(*track);
  }
}

