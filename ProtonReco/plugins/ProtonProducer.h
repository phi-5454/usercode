#ifndef ProtonProducer_H
#define ProtonProducer_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/CTPPSReco/interface/TotemRPCluster.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

#include <fstream>
#include <vector>

#include "../interface/gzstream.h"

namespace edm { class Event; class EventSetup; class Run; }

class RpReco;
class PrTrack;

class ProtonProducer : public edm::EDProducer
{
public:
  explicit ProtonProducer(const edm::ParameterSet& ps);
  ~ProtonProducer();
  virtual void produce(edm::Event& ev, const edm::EventSetup& es);

private:
  void beginJob();
  void endJob();

  void beginRun(edm::Run const & run, edm::EventSetup const & es);

  std::vector<PrTrack> processRpTracks(const edm::Event & ev);

  edm::EDGetTokenT<edm::DetSetVector<TotemRPCluster>> totemRpClusters_;
  edm::EDGetTokenT<edm::DetSetVector<TotemRPLocalTrack>> localTracks_;

  std::string outFile;

  const edm::ParameterSet pset;

  ogzstream file;

  RpReco * theRpReco;
};
#endif

