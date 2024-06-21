#ifndef _TStripHit_h_
#define _TStripHit_h_

#include <stdint.h>
#include <utility>
#include "TObject.h"

class TStripHit : public TObject
{
 public:
  TStripHit();
  virtual ~TStripHit();

  short int     meas;
  float         pred;

  bool forCalib; // isStraight && isComplete && hasBigPixelsOnlyInside
  bool forEloss; // isComplete

  float x;     // pathlength [cm]
  float thickness;

//  std::vector<uint8_t> adc;
  std::vector<unsigned char> adc;
  float previousGain;

  uint32_t detId;
  unsigned char chip;

  float r; // global r coordinate of the hit // TBR
  float p; // global momentum, updated state, at the hit // TBR

  ClassDef(TStripHit,1)
};

#endif
