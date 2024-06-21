#ifndef _SlimMeasurement_h_
#define _SlimMeasurement_h_

#include <stdint.h>
#include <utility>

class SlimMeasurement
{
 public:
  SlimMeasurement() {};
  virtual ~SlimMeasurement() {};

  unsigned short energy;
  unsigned short energySigma;

  unsigned short path;
  unsigned short epsilon;

//  unsigned short p;

  uint8_t det;

//  uint8_t nChannels;

  bool isSaturated;
};

#endif
