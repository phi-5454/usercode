#include "EnergyLossProducer.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"

#include "../interface/DataHandler.h"
#include "../interface/HadronClusters.h"
#include "../interface/TTrack.h"

using namespace std;

#undef DebugOutput

/*****************************************************************************/
EnergyLossProducer::EnergyLossProducer(const edm::ParameterSet& ps) : pset(ps)
{
  cerr << "\033[22;31m" << "Getting parameters.."
       << "\033[22;0m"  << endl;

  trackProducer = consumes<reco::TrackCollection>(
                  pset.getParameter<edm::InputTag>("trackProducer"));

  trajectoryProducer = consumes<vector<Trajectory> >(
                       ps.getParameter<edm::InputTag>("trajectoryProducer"));

  clusterShapes = consumes<SiPixelClusterShapeCache>(
                  ps.getParameter<edm::InputTag>("clusterShapes"));

  dedxProducer = consumes<reco::DeDxDataValueMap>(
                 pset.getParameter<edm::InputTag>("dedxProducer"));

  tag           = ps.getParameter<string>("tag");

  cerr << "\033[22;31m" << "Setting up products.."
       << "\033[22;0m"  << endl;
  produces<reco::DeDxDataValueMap>("energyLossPixHits");
  produces<reco::DeDxDataValueMap>("energyLossStrHits");
  produces<reco::DeDxDataValueMap>("energyLossAllHits");
}

/*****************************************************************************/
EnergyLossProducer::~EnergyLossProducer()
{
}

/*****************************************************************************/
void EnergyLossProducer::beginJob()
{
#ifdef DebugOutput
  file.open("result.dat");
#endif
}

/*****************************************************************************/
void EnergyLossProducer::endJob()
{
#ifdef DebugOutput
  file.close();
#endif
}

/*****************************************************************************/
void EnergyLossProducer::readFitResults()
{
  char fileName[256];
  sprintf(fileName,"UserCode/EnergyLossPID/data/mostprob_%s.dat", tag.c_str());

  edm::FileInPath fileInPath(fileName);
  ifstream file(fileInPath.fullPath().c_str());

  cerr << " from " << fileName << endl;

  while(!file.eof())
  {
    int icha,ieta,ipt, d;
    float mue, p;

    file >> icha >> ieta >> ipt >> p;

    for(int k = 0; k < K ; k++)
    {
      pair<pair<int,int>,pair<int,int> > key(pair<int,int>(icha,ieta),
                                             pair<int,int>( ipt,   k));

      file >> mean[key] >> mue >> amp[key];
    }

    pair<pair<int,int>,int> key(pair<int,int>(icha,ieta), ipt);

    file >> scale[key] >> d >> d;
  }

  file.close();
}

/*****************************************************************************/
void EnergyLossProducer::beginRun(edm::Run const & run, edm::EventSetup const & es)
{
  cerr << "\033[22;31m" << "Starting DataHandler (" << tag << ").."
       << "\033[22;0m"  << endl;

  theDataHandler = new DataHandler(tag);
  theDataHandler->applyGain = true;
  theDataHandler->beginJob();

  cerr << "\033[22;31m" << "Starting HadronClusters.."
       << "\033[22;0m"  << endl;
  theClusters = new HadronClusters(es, pset);

  // read fit results 
  cerr << "\033[22;31m" << "Reading fit amplitudes.."
       << "\033[22;0m"  << endl;
  readFitResults();
}

/*****************************************************************************/
pair<vector<float>, vector<float> > EnergyLossProducer::getProbabilities
   (const TTrack & track, const pair<double,double> & epsilon)
{
  pair<vector<float>,vector<float> > 
   res(vector<float>(3,0.), vector<float>(3,999.));

  if(track.eta > etaMin && track.eta < etaMax &&
     track.pt  >  ptMin && track.pt  <  ptMax)
  {  
    //
    int icha = (track.charge > 0 ? pos : neg);
    int ieta = int((track.eta - etaMin)/(etaMax - etaMin) * etaBins);
    int ipt  = int((track.pt  -  ptMin)/( ptMax -  ptMin) *  ptBins);

    //
    double    logde = log(epsilon.first);
    double siglogde = sqrt(epsilon.second) / epsilon.first;

    //
    pair<pair<int,int>,int> key(pair<int,int>(icha,ieta), ipt);

    double sig = scale[key] * siglogde;

    for(int k = 0; k < K; k++)
    {
      pair<pair<int,int>,pair<int,int> > key(pair<int,int>(icha,ieta),
                                             pair<int,int>( ipt,   k));
  
      double q = (logde - mean[key]) / sig;
  
      if(fabs(q) < 10) // FIXME
      {
        res.first[k]  = amp[key] * exp(-0.5*q*q); // division by sig not needed
        res.second[k] = q;
      }
    }

    double s = 0.;
      for(int k = 0; k < K; k++) s += res.first[k];

    if(s > 0.)
      for(int k = 0; k < K; k++) res.first[k] /= s;
  }

  return res;
}

/*****************************************************************************/
void EnergyLossProducer::produce(edm::Event& ev, const edm::EventSetup& es)
{
  //
  ev.getByToken(clusterShapes, clusterShapeCache);

  // Get track collection
  edm::Handle<reco::TrackCollection> trackHandle;
  ev.getByToken(trackProducer,       trackHandle);

  const vector<reco::Track> & trackCollection =
                                 *(trackHandle.product());

  // Get official dE/dx collection
  edm::Handle<reco::DeDxDataValueMap> dedxHandle;
  ev.getByToken(dedxProducer,         dedxHandle);

  auto outputPix = std::make_unique<reco::DeDxDataValueMap>(); 
  auto outputStr = std::make_unique<reco::DeDxDataValueMap>(); 
  auto outputAll = std::make_unique<reco::DeDxDataValueMap>(); 

  reco::DeDxDataValueMap::Filler fillerPix(*outputPix);
  reco::DeDxDataValueMap::Filler fillerStr(*outputStr);
  reco::DeDxDataValueMap::Filler fillerAll(*outputAll);

  // Get trajectory collection
  edm::Handle<vector<Trajectory> > trajeHandle;
  ev.getByToken(trajectoryProducer,     trajeHandle);
  const vector<Trajectory> & trajeCollection =
                           *(trajeHandle.product());

  vector<reco::DeDxData> estimatePix;
  vector<reco::DeDxData> estimateStr;
  vector<reco::DeDxData> estimateAll;

  // Take all trajectories
  int j = 0;
  for(vector<Trajectory>::const_iterator traje = trajeCollection.begin();
                                         traje!= trajeCollection.end();
                                         traje++, j++)
  {
    TTrack track;

    track.charge = trackCollection[j].charge();
    track.eta    = trackCollection[j].eta();
    track.pt     = trackCollection[j].pt();

    theClusters->analyzeRecTrajectory(*clusterShapeCache,*traje, track);

    pair<double,double> epsilon = 
      theDataHandler->processTrack(track, estimatePix,
                                          estimateStr,
                                          estimateAll);

    //
    pair<vector<float>, vector<float> > res = getProbabilities(track, epsilon);

    vector<float> pro = res.first;
    vector<float> dev = res.second;

#ifdef DebugOutput
    double p = trackCollection[j].pt() *
          cosh(trackCollection[j].eta());

    reco::TrackRef trackref = reco::TrackRef( trackHandle, j );

    const reco::DeDxData & dedxOfficial  = dedxHandle->get(trackref.key()); 

    // DataFormats/TrackReco/src/DeDxData.cc always returns dEdxError() == -1
    file << " " << p 
         << " " << estimateAll.back().dEdx()                 // epsilon
         << " " << sqrt(epsilon.second)                      // sigma(epsilon)
         << " " << estimateAll.back().numberOfMeasurements() // nhits
         << " " << pro[0] // prob pion
         << " " << pro[1] // prob kaon
         << " " << pro[2] // prob prot (sum = 1)
         << " " << dev[0] // deviation pion
         << " " << dev[1] // deviation kaon
         << " " << dev[2] // deviation prot (in units of sigma)
         << " " << dedxOfficial.dEdx()
         << " " << estimatePix.back().dEdx()                 // epsilon
         << " " << estimateStr.back().dEdx()                 // epsilon
         << " " << track.charge
         << endl; 
#endif
  }

  fillerPix.insert(trackHandle, estimatePix.begin(), estimatePix.end());
  fillerStr.insert(trackHandle, estimateStr.begin(), estimateStr.end());
  fillerAll.insert(trackHandle, estimateAll.begin(), estimateAll.end());

  fillerPix.fill();
  fillerStr.fill();
  fillerAll.fill();

  // Put back result to event
  ev.put(std::move(outputPix), "energyLossPixHits");
  ev.put(std::move(outputStr), "energyLossStrHits");
  ev.put(std::move(outputAll), "energyLossAllHits");
}

