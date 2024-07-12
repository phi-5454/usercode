#ifndef ProtonProducer_H
#define ProtonProducer_H

#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/CTPPSReco/interface/TotemRPCluster.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

#include <fstream>
#include <vector>

#include "../interface/gzstream.h"

namespace edm {
class Event;
class EventSetup;
class Run;
} // namespace edm

class RpReco;
class PrTrack;

class ProtonProducer : public edm::EDProducer {
public:
  explicit ProtonProducer(const edm::ParameterSet &ps);
  ~ProtonProducer();
  virtual void produce(edm::Event &ev, const edm::EventSetup &es);

private:
  void beginJob();
  void endJob();

  void beginRun(edm::Run const &run, edm::EventSetup const &es);

  std::vector<PrTrack> processRpTracks(const edm::Event &ev);

  edm::EDGetTokenT<edm::DetSetVector<TotemRPCluster>> totemRpClusters_;
  edm::EDGetTokenT<edm::DetSetVector<TotemRPLocalTrack>> localTracks_;

  std::string outFile;

  const edm::ParameterSet pset;

  ogzstream file;

  RpReco *theRpReco;

  // For output tree, output root file
  TString baseName;
  TString dirName;
  TFile *fOut;

  // PrTrack intermediate for the output
  std::vector<PrTrack> prTrks_out;

  // For run and event number. All else can be found from the other root files
  // (using a natural join on event id)
  UInt_t run_id;
  ULong64_t event_id;

  TTree *outT;
};
#endif
