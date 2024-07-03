#ifndef _RpFit_h_
#define _RpFit_h_

#include <map>
#include <vector>

#include "../interface/Structures.h" // FIXME needed for my_*

struct LocalFit;

//
class RpFit
{
 public:
  RpFit(const std::map<RpDet, std::map<std::vector<double>,int>> & patterns,
        const int nCmTra, int flag);

  ~RpFit(); 

  double getStrip(const RpDet & det, int pla, double a, double b);

  Vector2 getGlobalPosition(const RpDet & det,
                        const std::vector<Vector1> & lpos);

  LocalFit getLocalPosition(
    const RpDet & det, const std::vector<double> & strip,
    int ibas, int base, int step);

  std::vector<std::pair<double,double>> getPredicted(
    const RpDet & det, const std::vector<double> & strip);//, int move, int step);

  LocalFit fitTracklet(const RpDet & det, std::vector<double> strip);

  void optimizeShifts(bool collectChi2_);

  bool hasPattern(const RpDet & det, const std::vector<double> & strip);

  void fitPatterns();
  void readPatternFits();

 private:
  void getDet(const std::string & line, RpDet & det);

  void readShifts();

  bool isFilled(double y);
  bool isSingle(double y);

  double getTilt(const RpDet & det);
  Vector1 stripToPos(const Vector1 & strip_no);

  std::vector<double> getZ(const RpDet & det);

  double getDelta(double meas);

  double localChi2(const std::vector<double> & pars);

  double getMean(const std::vector<double> & strip);

  LocalFit fitThroughPolygon(const std::vector<double> & strip,
                             std::vector<std::vector<double>> & bands);
  std::vector<double> fitThroughSimplex(const std::vector<double> & strip);

  double globalChi2(const std::vector<double> & del);

  // patterns[det][{s0,s1,s2,s3,s4}] = num
  const std::map<RpDet, std::map<std::vector<double>,int>> & patterns;
  const int nCmTra;

  // shifts[det][pla]
  std::map<RpDet, std::vector<double>> shifts;

  // fits[det][strip] = LocalFit
  std::map<RpDet, std::map<std::vector<double>,LocalFit>> fits;

  //
  std::vector<double> my_strip, my_shift; // FIXME global
  RpDet my_det;                           // FIXME global

  //
  bool collectChi2; // collect chi2 or count zeros
  bool debug = false;

  bool doSquare = true;

  int iplot = 0;
};

#endif
