#ifndef _RpPrint_h_
#define _RpPrint_h_

#include <vector>
#include <string>

struct RpDet;
struct RpTrack;

//
class RpPrint
{
 public:
  void print(const std::vector<RpDet>   & rpDets,
             const std::vector<RpTrack> & rpTracks);

 private:
  std::string col(int col);
  void print(bool hasDet, bool hasTra, int uv);
};

#endif
