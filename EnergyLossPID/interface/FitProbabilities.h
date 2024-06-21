#ifndef _FitProbabilities_h_
#define _FitProbabilities_h_

#include <vector>
#include <map>

#include "TVectorD.h"
#include "TVirtualFitter.h"

class MostProbable;
class TVirtualFitter;

class FitProbabilities
{
 public:
  FitProbabilities(std::map<std::pair<int,std::pair<int,int> >,
                            std::vector<std::pair<double,double> > > & points_);

  virtual ~FitProbabilities();

  void run();

 private:
  void printPropeller(int i);

  double runLikelihood(int K,
    const std::vector<bool> & isFix,
    TVectorD & mu, TVectorD & P, TVectorD & mue);

  double getValue(double bg);

  void runLogLikelihood
    (const char* fileName,
     std::vector<double> & x,
     std::vector<double> & y,
     std::vector<double> & ye);

  void fitMostProbablePol3(const std::vector<double> & x,
                           const std::vector<double> & y,
                           const std::vector<double> & ye);

  MostProbable * mostProbable;
};

#endif
