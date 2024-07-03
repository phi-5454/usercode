#ifndef _RpEffic_h_
#define _RpEffic_h_

#include <vector>
#include <map>
#include <mutex>

//
struct RpTrack;
struct RpDet;

class RpPat;
class RpFit;
class Histo;

//
class RpEffic
{
 public:
  RpEffic(const int flag, const int nCmTra);
  ~RpEffic();

  bool collectStrip(const std::vector<RpTrack> & rpTracks, const int & run);

  double getStripEffic(int run, const RpDet & det);
  double getGroupEffic(int run, const RpDet & det, double Cy, double Cx);

 private:
  void writeStripEffic(const int & run);
  void writeStripEffic();

  void readStripEffic(const int & run);
  void readStripEffic();

  void readGroupEffic(const int & run);
  void readGroupEffic();

  // collectors (run,det)
  std::map<int, std::map<RpDet,double>> nFound, nEmpty;

  // strip efficiency vs (run,det)
  std::map<int, std::map<std::string,double>> stripEffic;

  // group efficiency vs (Cy,Cx)
  std::map<int, std::map<RpDet,Histo>> his_groupEffic;

  RpPat * theRpPat;
  RpFit * theRpFit;

  std::mutex mtx_hits, mtx_cout;

  //
  const int flag, nCmTra;
};

#endif
