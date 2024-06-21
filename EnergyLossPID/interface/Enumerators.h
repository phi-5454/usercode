#ifndef _Enumerators_h_
#define _Enumerators_h_

#include <cstring>
#include <vector>

enum AllDetector { PXB=0, PXF=1, TIB=2, TID=3, TOB=4, TEC3=5, TEC5=6, nDets };
static std::string detNames[nDets] =
                 {"PXB", "PXF", "TIB", "TID", "TOB", "TEC3", "TEC5"};

enum { pix=0, str=1, all=2, nVers=3 };
static std::string verNames[nVers] = {"pix", "str", "all"};

enum { pos=0, neg=1, nCharges=2 };
static std::string chargeNames[nCharges] = {"pos", "neg"};

enum { unknown=0, elec=1, pion=2, kaon=3, prot=4, nParts=5 };
static std::string partNames[nParts] = {"unknown","elec","pion","kaon","prot"};


#endif
