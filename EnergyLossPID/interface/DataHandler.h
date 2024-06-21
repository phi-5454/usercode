#ifndef _DataHandler_h_
#define _DataHandler_h_

#include <vector>
#include <utility>
#include <map>
#include <string>
#include <stdint.h>
#include <fstream>

#include "Enumerators.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

typedef std::pair<uint32_t, unsigned char> ChipId;

class     Measurement;
class SlimMeasurement;

class TBunchCrossing;
class TTrack;

class TFile;

class TH2S;
class TH2F;
class TH3I;

class MostProbable;
class ParticleType;
class EstimateHandler;

class DataHandler
{
 public:
  DataHandler();
  DataHandler(std::string tag_);
  virtual ~DataHandler();

  void options(int arg, char **arc, int chunk);

  void beginJob();
  void endJob(int filesFrom);

  int  getDetId(uint32_t id, float thickness);

  void processEvent (const TBunchCrossing * bunx);

  bool coupledHits, depositMap, calibGain, applyGain, surveyMap, fitProbs;

  int  pass;

  std::vector<TH2S> hcoupling;
  std::vector<TH3I> deposit3[nCharges];
  std::vector<TH3I> esigma3[nCharges];
  std::vector<TH2F> histos[nVers];

  void readGainCorrection();

  std::pair<double,double> processTrack(const TTrack & track,
                    std::vector<reco::DeDxData> & estimatePix,
                    std::vector<reco::DeDxData> & estimateStr,
                    std::vector<reco::DeDxData> & estimateAll);

  void fitProbabilities();

 private:
  void writeHistos(const std::vector<TH2S> & histos);
  void writeHistos(const std::vector<TH2F> & histos);
  void writeHistos(const std::vector<TH3I> & histos);

  void readStripProps();

  SlimMeasurement toSlim(const Measurement & meas); 
  void correctEnergy(float & energy, const ChipId & detId);
  double effectivePath(double l, int detType);  

  std::vector<Measurement> getPixelMeasurements(const TTrack & track);
  std::vector<Measurement> getStripMeasurements(const TTrack & track);
                    void getCoupledMeasurements(const TTrack & track);

  void fillHits(std::vector<Measurement> & meas,
                int cha, double bg);

  void copyEpsilon(std::vector<Measurement> & meas, int pid);

  void fillHistos(const TTrack & track,
                  const std::vector<int> & nhits,
                  const std::vector<std::pair<double,double> > & values,
                  const std::vector<Measurement> & measAll);

  void processTrack(const TTrack & track);

  MostProbable * mostProbable;
  ParticleType * particleType;
  EstimateHandler * estimateHandler;

  std::vector<double> thr;
  std::vector<double> alpha;
  std::vector<double> sigma;

  std::map<ChipId  , float> gain;
  std::map<ChipId, std::vector<SlimMeasurement> > hits;

  int nrec;

  std::map<uint32_t, float> thickness;

  std::vector<std::vector<double> > epsilon;

  std::map<std::pair<int,std::pair<int,int> >,
           std::vector<std::pair<double,double> > > points;

  std::string tag;
};

#endif
