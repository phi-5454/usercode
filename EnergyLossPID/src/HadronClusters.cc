#include "../interface/HadronClusters.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/Common/interface/Handle.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"

#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShape.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterData.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"

#include "RecoLocalTracker/SiPixelClusterizer/plugins/PixelThresholdClusterizer.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationOfflineService.h"

#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelClusterShapeCache.h"

#include "../interface/TTrack.h"

#include <utility>
#include <vector>
#include <fstream>

using namespace std;

/*****************************************************************************/
HadronClusters::HadronClusters(const edm::EventSetup& es,
                               const edm::ParameterSet& pset)
{
  // Get tracker geometry
  edm::ESHandle<TrackerGeometry>          tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);
  theTracker =                            tracker.product();

////  theClusterShape = new ClusterShapeHitFilter(es);
  edm::ESHandle<ClusterShapeHitFilter> shape;
  es.get<CkfComponentsRecord>().get("ClusterShapeHitFilter",shape);
  theClusterShape = shape.product();

  theSiPixelGainCalibration_ = new SiPixelGainCalibrationOfflineService(pset);
  setSiPixelGainCalibrationService(theSiPixelGainCalibration_);
  theSiPixelGainCalibration_->setESObjects( es );
}

/*****************************************************************************/
HadronClusters::~HadronClusters()
{
}

/*****************************************************************************/
// Strip
void HadronClusters::processRec(const SiStripRecHit2D & recHit,
     LocalVector ldir, TTrack & r, float gmom, float tpos)
{
  TStripHit hit;

  int meas;

  // Complete
  LocalPoint lpos;
  hit.forCalib = theClusterShape->getSizes(recHit,lpos,ldir, meas,hit.pred);

  // Complete and compatible
  hit.forEloss = ( hit.forCalib &&
                   theClusterShape->isCompatible(recHit,ldir) );

  hit.meas = meas;

  DetId id = recHit.geographicalId();
  const StripGeomDetUnit* stripDet =
    dynamic_cast<const StripGeomDetUnit*> (theTracker->idToDet(id));

  // Collect adc
  for(vector<uint8_t>::const_iterator
    i = (recHit.cluster()->amplitudes()).begin();
    i!= (recHit.cluster()->amplitudes()).end(); i++)
  {
    hit.adc.push_back(*i);
  }

  // Length
  hit.x = ldir.mag()/fabsf(ldir.z()) *
          stripDet->surface().bounds().thickness();
  hit.thickness =
          stripDet->surface().bounds().thickness();

  hit.detId = recHit.geographicalId();

  int ix = int(recHit.cluster()->barycenter() / 128.);
  hit.chip = ix;

  hit.r = tpos;
  hit.p = gmom;

  r.stripHits.push_back(hit);
}

/*****************************************************************************/
// Pixel
void HadronClusters::processRec
    (const SiPixelClusterShapeCache & clusterShapeCache,
     const SiPixelRecHit & recHit,
     LocalVector ldir, TTrack & r, float gmom, float tpos)
{
  //
  TPixelHit hit;

  int part;
//  vector<pair<int,int> > meas;
  ClusterData::ArrayType meas;

  theClusterShape->getSizes(recHit,ldir, clusterShapeCache, part,meas, hit.pred);

  if(meas.size() > 0)
    hit.meas = meas[0]; // FIXME

  DetId id = recHit.geographicalId();
  const PixelGeomDetUnit* pixelDet =
    dynamic_cast<const PixelGeomDetUnit*> (theTracker->idToDet(id));

  hit.isSaturated = false;

  float MeVperElec = 3.61e-6;

  double delta = 0.;

  // Collect adc
  int ii = 0;
  for(vector<uint16_t>::const_iterator
        i = (recHit.cluster()->pixelADC()).begin();
        i!= (recHit.cluster()->pixelADC()).end(); i++)
  {
    // ( pos , (x , y) )
    uint8_t x = recHit.cluster()->pixelOffset()[ii++];
    uint8_t y = recHit.cluster()->pixelOffset()[ii++];

    DetId detid_ = id;
    int row = x + recHit.cluster()->minPixelRow();
    int col = y + recHit.cluster()->minPixelCol();

    // Go back to adc 
    float DBgain     =
      theSiPixelGainCalibrationService_->getGain(detid_, col, row);
    float DBpedestal = 
      theSiPixelGainCalibrationService_->getPedestal(detid_, col, row) * DBgain;

    if(DBgain > 0.)
    {
      int    elec = int(*i);
      double vcal = (elec + 414.)/65.; // was 65.5
      double adc  = (vcal + DBpedestal)/DBgain;

/*
cerr << " pix " << *i << " " << elec << " " << vcal << " " << adc
     << " : "  << DBgain << " " << DBpedestal
     << " " << int(id) << " " << theTracker->idToDet(id)->subDetector()
     << endl;
*/

      if(adc > 185) hit.isSaturated = true; // was 200
    }

    double de = *i * MeVperElec;
    delta += de;
  }

  hit.Delta = delta;

  // Length
  hit.x = ldir.mag()/fabsf(ldir.z()) *
          pixelDet->surface().bounds().thickness();
  hit.thickness =
          pixelDet->surface().bounds().thickness();

  hit.detId = recHit.geographicalId();

  int ix = int(recHit.cluster()->x() / 80.);
  int iy = int(recHit.cluster()->y() / 52.);
  hit.chip = (ix << 3) + iy;

  ClusterData data;
  ClusterShape theClusterShapeEE;
  theClusterShapeEE.determineShape(*pixelDet, recHit, data);

  // Complete, straight and no big pixels outside
  hit.forCalib = ( data.isComplete &&
                   data.isStraight &&
                   data.hasBigPixelsOnlyInside );

  // Complete and compatible
  hit.forEloss = ( data.isComplete &&
                   theClusterShape->isCompatible(recHit,ldir, clusterShapeCache) );

  hit.r = tpos;
  hit.p = gmom;

  hit.nChannels = (recHit.cluster()->pixelADC()).size();

  r.pixelHits.push_back(hit);
}

/*****************************************************************************/
void HadronClusters::analyzeRecTrajectory
  (const SiPixelClusterShapeCache & clusterShapeCache,
   const Trajectory & trajectory,
   TTrack & r)
{
  for(vector<TrajectoryMeasurement>::const_iterator
      meas = trajectory.measurements().begin();
      meas!= trajectory.measurements().end(); meas++)
  {
    const TrackingRecHit* recHit = meas->recHit()->hit();
    DetId id = recHit->geographicalId();

    const StripGeomDetUnit* stripDet =
     dynamic_cast<const StripGeomDetUnit*> (theTracker->idToDet(id));

    if(recHit->isValid())
    {
      LocalVector ldir = meas->updatedState().localDirection();
      float       gmom = meas->updatedState().globalMomentum().mag();
      float       tpos = meas->updatedState().globalPosition().perp();

      if(theTracker->idToDet(id)->subDetector() == GeomDetEnumerators::P1PXB ||
         theTracker->idToDet(id)->subDetector() == GeomDetEnumerators::P1PXEC)
      {
        // Pixel
        const SiPixelRecHit* pixelRecHit =
          dynamic_cast<const SiPixelRecHit *>(recHit);

        if(pixelRecHit != 0)
          processRec(clusterShapeCache, *pixelRecHit, ldir, r, gmom,tpos);
      }
      else
      {
        // Strip
        const SiStripRecHit1D* stripSimpleHit  =
          dynamic_cast<const SiStripRecHit1D*>(recHit);

        const SiStripMatchedRecHit2D* stripMatchedRecHit =
          dynamic_cast<const SiStripMatchedRecHit2D *>(recHit);

        const ProjectedSiStripRecHit2D* stripProjectedRecHit =
          dynamic_cast<const ProjectedSiStripRecHit2D *>(recHit);

        const SiStripRecHit2D* stripRecHit = 
          dynamic_cast<const SiStripRecHit2D *>(recHit);

        if(stripSimpleHit != 0)
        {
          SiStripRecHit2D stripHit(
                               stripSimpleHit->localPosition(),
                    LocalError(stripSimpleHit->localPositionError().xx(),0.,
                        std::numeric_limits<float>::max()),
//                        stripSimpleHit->geographicalId(), 
//theTracker->idToDet(id),
*stripDet,
                        stripSimpleHit->cluster()
                    );

          processRec(stripHit, ldir, r, gmom,tpos);
        }

        if(stripMatchedRecHit != 0)
        {
          processRec((stripMatchedRecHit->monoHit())  , ldir, r, gmom,tpos);
          processRec((stripMatchedRecHit->stereoHit()), ldir, r, gmom,tpos);
        }

        if(stripProjectedRecHit != 0)
          processRec(stripProjectedRecHit->originalHit(), ldir, r, gmom,tpos);

        if(stripRecHit != 0)
          processRec(*stripRecHit, ldir, r, gmom,tpos);
      }
    }
  }
}

