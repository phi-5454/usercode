#include "ProtonProducer.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

//
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"
#include "DataFormats/CTPPSReco/interface/TotemRPCluster.h"

//
#include "../interface/Structures.h"
#include "../interface/RpReco.h"

//
using namespace std;

#define DebugOutput

/*****************************************************************************/
ProtonProducer::ProtonProducer(const edm::ParameterSet& ps) : pset(ps)
{
  cerr << "\033[22;31m" << "Getting parameters.."
       << "\033[22;0m"  << endl;

  totemRpClusters_ = consumes<edm::DetSetVector<TotemRPCluster>>(
                     ps.getParameter<edm::InputTag>("totemRpClusters"));

  localTracks_  = consumes<edm::DetSetVector<TotemRPLocalTrack>>(
                  ps.getParameter<edm::InputTag>("localTracks"));

  outFile = ps.getParameter<string>("outFile");

  //
  theRpReco = new RpReco(2, true); // nCmTra, read RpAlign
}

/*****************************************************************************/
ProtonProducer::~ProtonProducer()
{
}

/*****************************************************************************/
void ProtonProducer::beginJob()
{
#ifdef DebugOutput
  file.open(outFile.c_str());
#endif
}

/*****************************************************************************/
void ProtonProducer::endJob()
{
#ifdef DebugOutput
  file.close();
#endif
}

/*****************************************************************************/
void ProtonProducer::beginRun(edm::Run const & run, edm::EventSetup const & es)
{
}

/*****************************************************************************/
vector<PrTrack> ProtonProducer::processRpTracks(const edm::Event & ev)
{
  edm::Handle<edm::DetSetVector<TotemRPCluster> > stripClusters;
  ev.getByToken(                totemRpClusters_, stripClusters);

  edm::Handle<edm::DetSetVector<TotemRPLocalTrack>> localTracks;
  ev.getByToken(                      localTracks_, localTracks);

  //
/*
  for(auto & tracks : *localTracks)
  {
    TotemRPDetId detId(tracks.detId());

    if(detId.station() == 0 ||
       detId.station() == 2)
    {
      file << " " << detId.arm()
           << " " << detId.station()
           << " " << detId.rp();

      for(auto & track : tracks)
      {
        vector<double> center(10,-1);

        for(auto & recHits : track.getHits())
        for(auto & recHit : recHits)
        {
          TotemRPDetId detId(recHits.detId());
          const TotemRPCluster & cluster = recHit.getCluster();

          center[detId.plane()] = cluster.getCenterStripPosition();
        }

        for(auto & c : center)
          file << " " << c;
        file << endl;
      }
    }
  }
*/

  //
  Event event;

  // roman pot tracks
  for(auto & tracks : *localTracks)
  {
    TotemRPDetId detId(tracks.detId());

    if(detId.station() == 0 ||
       detId.station() == 2)
    {
      RpTrack rpTrack;
      RpDet   & det = rpTrack.det;
//    Vector2 & pos = rpTrack.pos;

      det.arm = detId.arm();
      det.sta = detId.station();
      det.rpt = detId.rp();

      bool ok = (det.sta == 0 || det.sta == 2) ||
                (det.rpt == 4 || det.rpt == 5);

      det.sta /= 2; // 0 2 -> 0 1
      det.rpt %= 4; // 4 5 -> 0 1

      for(auto & track : tracks)
      {
        vector<double> center(10,-1);

        for(auto & recHits : track.getHits())
        for(auto & recHit : recHits)
        {
          TotemRPDetId detId(recHits.detId());
          const TotemRPCluster & cluster = recHit.getCluster();

          center[detId.plane()] = cluster.getCenterStripPosition();
        }

        int c = 0;
        for(short int pla = 0; pla < nPlanes; pla++)
        for(short int uv = 0; uv < 2; uv++)
        {
          double strip = center[c++];

          if(det.print() == "0|0|0" && uv == 0 && pla == 2)
            strip = -1;

          if(strip == -1) strip = Empty; // -1 -> -99

          rpTrack.clus[uv][pla] = strip;
        }
      }

      //
      if(ok) event.rpTracks.push_back(rpTrack);
    }
  }

  //
  vector<PrTrack> prTracks;

  // need exactly four tracklets, proper topology
  if(theRpReco->getTopo(event.rpTracks, event.topo))
  {
    theRpReco->reconstruct(event);

    prTracks = event.prTracks;
  } // otherwise prTracks is left empty

  return prTracks;
}

/*****************************************************************************/
void ProtonProducer::produce(edm::Event& ev, const edm::EventSetup& es)
{
  vector<PrTrack> prTracks = processRpTracks(ev);

  // FOR NOW, WRITING OUT RECONSTRUCTED PROTON KINEMATICS TO out.dat.gz

#ifdef DebugOutput
  if(prTracks.size() > 0)
  { // write out protons
    for(int i = 0; i < 2; i ++)
      file << " " << prTracks[i].p.x // momentum components
           << " " << prTracks[i].p.y
           << " " << prTracks[i].p.z
           << " " << prTracks[i].pt.cxx // covariances of px and py
           << " " << prTracks[i].pt.cxy
           << " " << prTracks[i].pt.cyy
           << " " << prTracks[i].pos.x  // transverse location at IP
           << " " << prTracks[i].pos.y;

    file << endl;
  }
  else
    file << " -1" << endl;
#endif
    
}

