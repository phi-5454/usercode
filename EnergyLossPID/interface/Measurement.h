#ifndef _Measurement_h_
#define _Measurement_h_

#include <stdint.h>
#include <utility>

typedef std::pair<uint32_t, unsigned char> ChipId;

class Measurement
{
 public:
  Measurement() {};
  virtual ~Measurement() {};

  float energy;
  float energySigma;

  float path;   // effective path, l*(1 + a*log(l/l0)
//  float log_path_l0;
  float epsilon;

  float r;
  float p;

  bool isSaturated; 

  ChipId chipId;
  uint8_t det;

  uint8_t nChannels;
};

#endif
