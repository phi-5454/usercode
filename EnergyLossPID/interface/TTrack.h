#ifndef _TTrack_h_
#define _TTrack_h_

#include <utility>
#include "TObject.h"

#include "TPixelHit.h"
#include "TStripHit.h"

class TTrack : public TObject
{
 public:
  TTrack();
  virtual ~TTrack();

  short int charge;

  float eta;
  float pt;
  float phi;

  float chi2;
  short int ndf;

//  float z;
//  float d0;

  std::vector<TPixelHit> pixelHits;
  std::vector<TStripHit> stripHits;

  bool /*isPrimary,*/ isHighPurity;

//  short int algo;

  ClassDef(TTrack,1)
};

#endif
