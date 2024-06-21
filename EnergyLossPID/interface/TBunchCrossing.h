#ifndef _TBunchCrossing_h_
#define _TBunchCrossing_h_

#include "TObject.h"
#include <vector>
#include <string>

class TVertex;

class TBunchCrossing : public TObject
{
 public:
  TBunchCrossing();
  virtual ~TBunchCrossing();

  int runNumber;
  int lumiSection;
  int bxNumber;

  // HF
//  std::pair<short int,short int> hfTowers;

  // beam spot
//  double beamSpot_pos[3];
//  double beamSpot_sig[3];

  // vertices
  std::vector<TVertex> recVertices;

  void Clear();

  ClassDef(TBunchCrossing,1)
};

#endif
