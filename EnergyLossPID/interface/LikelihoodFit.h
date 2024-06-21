#ifndef _LikelihoodFit_h_
#define _LikelihoodFit_h_

#include <vector>
#include <utility>

#include <cmath>

class Measurement;
namespace CLHEP { class HepVector; class HepMatrix; }

class LikelihoodFit
{
 public:
  LikelihoodFit(const double & sigma0_,
                const double & b_) : sigma0(sigma0_), b(b_)
  {
    // Define constants
    nu     = 0.65;
    sqrtnu = sqrt(nu);
  }

  double probability(const Measurement & m, double epsilon);

  void calculate_wrt_y
    (const Measurement & m,                 double * val);
  void calculate_wrt_epsilon
    (const Measurement & m, double epsilon, double * val);

 private:
  double sigmaD(double y);

  double nu, c, sqrtnu;

  const double sigma0, b;
};

#endif
