#ifndef _RpPat_h_
#define _RpPat_h_

#include <vector>
#include <map>
#include <mutex>
#include <cmath>

struct RpTrack;
struct RpDet;

class RpFit;

//
class RpPat
{
 const double minFraction = 1e-3; // FIXME VERY

 public:
  RpPat(const int flag, const int nCmTra);
  ~RpPat();

  const std::map<RpDet, std::map<std::vector<double>,int>> & getPatterns(
    bool addHoles);

  void collectPatterns(const std::vector<RpTrack> & rpTracks);

  std::vector<int> toPattern(std::vector<double> & strip);

  void collectForEffic(RpDet det, const std::vector<double> & str,
                       RpFit * theRpFit,
                       std::map<RpDet,double> & nFound,
                       std::map<RpDet,double> & nEmpty);

 private:
  void writePatterns();
  void writeOccupancy();

  void getDet(const std::string & line, RpDet & det);

  void readPatterns(int turn);

  bool isFilled(double x);
  int nFilled(const double strip[]);

  bool isSingle(double x);

  std::vector<double> rebase(int ibas, const double clus[]);
  std::vector<double> rebase(int ibas, const std::vector<double> & clus);

  int removeShear(int ibas, std::vector<double> & str);

  int getBase(const std::vector<double> & str);

  void getSubsets(const std::vector<double> & str,
              std::vector<std::vector<int>> & sub);

  int reduce(std::vector<double> & str);

  // toCol
  std::map<RpDet, std::map<std::vector<double>,int>> patterns;
  std::map<RpDet, std::map<int, std::map<double,int>>> occupancy;

  std::mutex mtx_pats , mtx_occup;
  std::mutex mtx_found, mtx_empty;

  std::map<RpDet,int> highest;

  int sFound, sEmpty;

  const int flag, nCmTra;
};

#endif
