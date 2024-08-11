#include "ProtonProducer.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/FWLite/interface/Handle.h"

//
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/TotemRPCluster.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"

// For exprting products. TODO: Replace with more complete descriptor
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"

//
#include "../interface/RpReco.h"
#include "../interface/Structures.h"

//
using namespace std;

#define DebugOutput

/*****************************************************************************/
ProtonProducer::ProtonProducer(const edm::ParameterSet &ps) : pset(ps) {
  cerr << "\033[22;31m" << "Getting parameters.." << "\033[22;0m" << endl;

  totemRpClusters_ = consumes<edm::DetSetVector<TotemRPCluster>>(
      ps.getParameter<edm::InputTag>("totemRpClusters"));

  localTracks_ = consumes<edm::DetSetVector<TotemRPLocalTrack>>(
      ps.getParameter<edm::InputTag>("localTracks"));

  // produces<vector<RpTrack>>("one");
  outFile = ps.getParameter<string>("outFile");

  //
  theRpReco = new RpReco(2, true); // nCmTra, read RpAlign
}

/*****************************************************************************/
ProtonProducer::~ProtonProducer() {}

/*****************************************************************************/
void ProtonProducer::beginJob() {
  baseName = gSystem->BaseName(outFile.c_str());
  dirName = gSystem->DirName(outFile.c_str());
  // fOut = TFile::Open(dirName + "/" + baseName, "RECREATE");
  fOut = TFile::Open(outFile.c_str(), "RECREATE");
  fOut->cd();

  prTrks_out = vector<PrTrack>(2);

  // BOOK OUTPUT TREE
  outT = new TTree("tree", "tree");

  /*
      file << " " << prTracks[i].p.x // momentum components
           << " " << prTracks[i].p.y << " " << prTracks[i].p.z << " "
           << prTracks[i].pt.cxx // covariances of px and py
           << " " << prTracks[i].pt.cxy << " " << prTracks[i].pt.cyy << " "
           << prTracks[i].pos.x // transverse location at IP
           << " " << prTracks[i].pos.y;
   */

  // Bind to tree variables
  // Event info
  outT->Branch("Run", &run_id, "Run/i");
  outT->Branch("EventNum", &event_id, "EventNum/l");

  // (IP) Momentum
  outT->Branch("pr_px_a", &prTrks_out[0].p.x, "pr_px_a/D");
  outT->Branch("pr_px_b", &prTrks_out[1].p.x, "pr_px_b/D");

  outT->Branch("pr_py_a", &prTrks_out[0].p.y, "pr_py_a/D");
  outT->Branch("pr_py_b", &prTrks_out[1].p.y, "pr_py_b/D");

  outT->Branch("pr_pz_a", &prTrks_out[0].p.z, "pr_pz_a/D");
  outT->Branch("pr_pz_b", &prTrks_out[1].p.z, "pr_pz_b/D");
  // No covariance values found for Vector3

  // (IP) Transverse momentum
  outT->Branch("pr_ptx_a", &prTrks_out[0].pt.x, "pr_ptx_a/D");
  outT->Branch("pr_ptx_b", &prTrks_out[1].pt.x, "pr_ptx_b/D");

  outT->Branch("pr_pty_a", &prTrks_out[0].pt.y, "pr_pty_a/D");
  outT->Branch("pr_pty_b", &prTrks_out[1].pt.y, "pr_pty_b/D");

  // Errors in transverse momentum
  outT->Branch("pr_ptx_sigma_a", &prTrks_out[0].pt.cxx, "pr_ptx_sigma_a/D");
  outT->Branch("pr_ptx_sigma_b", &prTrks_out[1].pt.cxx, "pr_ptx_sigma_b/D");

  outT->Branch("pr_pty_sigma_a", &prTrks_out[0].pt.cyy, "pr_pty_sigma_a/D");
  outT->Branch("pr_pty_sigma_b", &prTrks_out[1].pt.cyy, "pr_pty_sigma_b/D");

  // (IP) Positions
  outT->Branch("pr_posx_a", &prTrks_out[0].pos.x, "pr_posx_a/D");
  outT->Branch("pr_posx_b", &prTrks_out[1].pos.x, "pr_posx_b/D");

  outT->Branch("pr_posy_a", &prTrks_out[0].pos.y, "pr_posy_a/D");
  outT->Branch("pr_posy_b", &prTrks_out[1].pos.y, "pr_posy_b/D");

  // Errors in position (IP)
  outT->Branch("pr_posx_sigma_a", &prTrks_out[0].pos.cxx, "pr_posx_sigma_a/D");
  outT->Branch("pr_posx_sigma_b", &prTrks_out[1].pos.cxx, "pr_posx_sigma_b/D");

  outT->Branch("pr_posy_sigma_a", &prTrks_out[0].pos.cyy, "pr_posy_sigma_a/D");
  outT->Branch("pr_posy_sigma_b", &prTrks_out[1].pos.cyy, "pr_posy_sigma_b/D");

  //

#ifdef DebugOutput
  // file.open(outFile.c_str());
#endif
}

/*****************************************************************************/
void ProtonProducer::endJob() {
  // save histos to file
  std::cout << endl
            << "Writes " << fOut->GetName() << " with " << outT->GetEntries()
            << " events." << endl;
  fOut->cd();
  outT->Write();

  fOut->Close();
#ifdef DebugOutput
  // file.close();
#endif
}

/*****************************************************************************/
void ProtonProducer::beginRun(edm::Run const &run, edm::EventSetup const &es) {}

/*****************************************************************************/
vector<PrTrack> ProtonProducer::processRpTracks(const edm::Event &ev) {
  edm::Handle<edm::DetSetVector<TotemRPCluster>> stripClusters;
  ev.getByToken(totemRpClusters_, stripClusters);

  edm::Handle<edm::DetSetVector<TotemRPLocalTrack>> localTracks;
  ev.getByToken(localTracks_, localTracks);

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
  for (auto &tracks : *localTracks) {
    TotemRPDetId detId(tracks.detId());

    if (detId.station() == 0 || detId.station() == 2) {
      RpTrack rpTrack;
      RpDet &det = rpTrack.det;
      //    Vector2 & pos = rpTrack.pos;

      det.arm = detId.arm();
      det.sta = detId.station();
      det.rpt = detId.rp();

      bool ok =
          (det.sta == 0 || det.sta == 2) || (det.rpt == 4 || det.rpt == 5);

      det.sta /= 2; // 0 2 -> 0 1
      det.rpt %= 4; // 4 5 -> 0 1

      for (auto &track : tracks) {
        vector<double> center(10, -1);

        for (auto &recHits : track.getHits())
          for (auto &recHit : recHits) {
            TotemRPDetId detId(recHits.detId());
            const TotemRPCluster &cluster = recHit.getCluster();

            center[detId.plane()] = cluster.getCenterStripPosition();
          }

        int c = 0;
        for (short int pla = 0; pla < nPlanes; pla++)
          for (short int uv = 0; uv < 2; uv++) {
            double strip = center[c++];

            if (det.print() == "0|0|0" && uv == 0 && pla == 2)
              strip = -1;

            if (strip == -1)
              strip = Empty; // -1 -> -99

            rpTrack.clus[uv][pla] = strip;
          }
      }

      //
      if (ok)
        event.rpTracks.push_back(rpTrack);
    }
  }

  //
  vector<PrTrack> prTracks;

  // need exactly four tracklets, proper topology
  if (theRpReco->getTopo(event.rpTracks, event.topo)) {
    theRpReco->reconstruct(event);

    prTracks = event.prTracks;
  } // otherwise prTracks is left empty

  return prTracks;
}

/*****************************************************************************/
void ProtonProducer::produce(edm::Event &ev, const edm::EventSetup &es) {
  vector<PrTrack> prTracks = processRpTracks(ev);
  // auto op = CTPPSLocalTrackLite(prTrackL.);

  auto outputPrs = std::make_unique<std::vector<PrTrack>>();
  for (auto &&pr : prTracks) {

    // auto pt = sqrt(pr.pt.x * pr.pt.x + pr.pt.y * pr.pt.y);
    //  TODO: Rigor. Uncertainty likely doesn't transform as euclidean distance
    //  does
    // auto ptu = sqrt(pr.pt.cxx * pr.pt.cxx + pr.pt.cyy * pr.pt.cyy);

    // Use covariances for uncertainties
    // outputPrs->push_back(CTPPSLocalTrackLite(-1, pr.pos.x, pr.pos.cxx,
    // pr.pos.y, pr.pos.cyy, pt, ptu));
    outputPrs->push_back(pr);
  }

  // write out protons
  if (prTracks.size() > 0) {
    event_id = ev.id().event();
    run_id = ev.id().run();

    // Try to make sure no memory is moved, only copied
    // prTrks_out = prTrks;
    for (int i = 0; i < 2; i++) {
      prTrks_out[i].p = prTracks[i].p;
      prTrks_out[i].pt = prTracks[i].pt;
      prTrks_out[i].pos = prTracks[i].pos;
    }

    // Save output tree
    outT->Fill();
  }

  /*
  std::cerr << "AAA\n";
  if (prTracks.size() > 0) {
    for (int i = 0; i < 2; i++) {
      std::cerr << " " << prTracks[i].p.x // momentum components
                << " " << prTracks[i].p.y << " " << prTracks[i].p.z << " "
                << prTracks[i].pt.cxx // covariances of px and py
                << " " << prTracks[i].pt.cxy << " " << prTracks[i].pt.cyy << " "
                << prTracks[i].pos.x // transverse location at IP
                << " " << prTracks[i].pos.y;
      std::cerr << endl;
    }
  }
  */

  // FOR NOW, WRITING OUT RECONSTRUCTED PROTON KINEMATICS TO out.dat.gz
  // std::cerr << "SMTH";

  /*
#ifdef DebugOutput
      if (prTracks.size() > 0) { // write out protons
    for (int i = 0; i < 2; i++) {
      file << " " << prTracks[i].p.x // momentum components
           << " " << prTracks[i].p.y << " " << prTracks[i].p.z << " "
           << prTracks[i].pt.cxx // covariances of px and py
           << " " << prTracks[i].pt.cxy << " " << prTracks[i].pt.cyy << " "
           << prTracks[i].pos.x // transverse location at IP
           << " " << prTracks[i].pos.y;
      file << endl;
      std::cerr << " " << prTracks[i].p.x // momentum components
                << " " << prTracks[i].p.y << " " << prTracks[i].p.z << " "
                << prTracks[i].pt.cxx // covariances of px and py
                << " " << prTracks[i].pt.cxy << " " << prTracks[i].pt.cyy <<
" "
                << prTracks[i].pos.x // transverse location at IP
                << " " << prTracks[i].pos.y;
      std::cerr << endl;
    }
  }
  else file << " -1" << endl;
#endif
  */

  // Add our new info to the event
  // ev.put(std::move(outputPrs), "one");
}
