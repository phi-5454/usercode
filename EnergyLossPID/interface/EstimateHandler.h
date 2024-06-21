#ifndef _EstimateHandler_h_
#define _EstimateHandler_h_

#include <iostream>
#include <stdint.h>
#include <vector>
#include <map>

class Measurement;
class SlimMeasurement;

typedef std::pair<uint32_t, unsigned char> ChipId;

class EstimateHandler
{
 public:
  EstimateHandler();

  void calibrateGain
    (std::map<ChipId, std::vector<SlimMeasurement> > & hits,
     std::map<uint32_t, float> thickness,
     char * fileName,
     std::map<ChipId, float> & gain);

  // basic estimator
  double estimate_(std::vector<Measurement> & meas,
                   std::pair<double,double> & value);

  // with removal
  double estimate (std::vector<Measurement> & meas,
                   std::pair<double,double> & value);

 private:
  int getDetId(uint32_t id, float thickness);

  bool allSaturated(const std::vector<Measurement> & meas);

  double sigma0;
  double b;
};

#endif
