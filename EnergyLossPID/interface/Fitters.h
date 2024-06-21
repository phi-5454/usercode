#ifndef _Fitters_h_
#define _Fitters_h_

#include <vector>
#include <utility>

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"

class     Measurement;
class SlimMeasurement;

class Fitters
{
 public:
  Fitters(const double sigma0_,
          const double b_,
          const std::vector<SlimMeasurement> & slim_);
  Fitters(const double sigma0_,
          const double b_,
          const std::vector<Measurement>     & meas_);

  virtual ~Fitters();

  double scanModule  (/*double function(double) */);
  double goldenSearch(/*double function(double),*/
                      double ax, double bx, double cx);

  double newtonMethodGain               (std::pair<double,double>   & value);
  double newtonMethodGain2  (std::vector<std::pair<double,double> > & values);
  
  double minimizeAllSaturated(std::pair<double,double> & value, int & nStep);
  double newtonMethodEpsilon (std::pair<double,double> & value, int & nStep);

 private:
  Measurement toMeas(const SlimMeasurement & slim);

  void shft2(double &a, double &b, const double c);
  void shft3(double &a, double &b, double &c, const double d);

  void   functionGain (double gain   , std::vector<double> & val);
  double functionGain2(const std::vector<double> & x,
                       CLHEP::HepMatrix & alpha,
                       CLHEP::HepVector & beta);
  void   functionEpsilon(double epsilon, std::vector<double> & val);

  double functionGain(double gain);

  const double sigma0,b;

  const std::vector<SlimMeasurement> slim;
  const std::vector<Measurement>     meas;
};

#endif
