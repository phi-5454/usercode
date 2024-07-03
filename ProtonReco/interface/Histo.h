#ifndef _Histo_h_
#define _Histo_h_

#include <fstream>
#include <string>
#include <vector>
#include <mutex>

class Random;
class igzstream;
class ogzstream;

#ifndef _Axis_
#define _Axis_
struct Axis {
  float min,max;
  int bins;

  float center(int ix) const
  { return min + (max - min)/bins * (ix+0.5); }

  float binwidth() const
  { return (max - min)/bins; }
};
#endif

class Histo
{
 public:
  Histo();
  ~Histo(); 

  // axes
  std::vector<Axis> axes;

  // init
  void init(float xmin, float xmax, int xbins,
            std::string name, bool norm = false, bool zip = false);

  void init(float xmin, float xmax, int xbins,
            float ymin, float ymax, int ybins,
            std::string name, bool norm = false, bool zip = false);

  void init(float xmin, float xmax, int xbins,
            float ymin, float ymax, int ybins,
            float zmin, float zmax, int zbins,
            std::string name, bool norm = false, bool zip = false);

  void init(float xmin, float xmax, int xbins,
            float ymin, float ymax, int ybins,
            float zmin, float zmax, int zbins,
            float vmin, float vmax, int vbins,
            std::string name, bool norm = false, bool zip = false,
                                                 bool uchar = false);

  // get
  std::pair<float,float> & get(const std::vector<int>  & ix);
  float                    val(const std::vector<float> & x);
  std::pair<float,float>  getv(const std::vector<float> & x);

  // fill
  void fillw(const std::vector<float> & x, float w);
  void fill (const std::vector<float> & x);

  // set value
  void set  (const std::vector<float> & x, float v);

  // divide by another Histogram
  void div(const Histo & q);

  // write
  void write();

  // read
  void read();

  // sample
  float sample(const std::vector<float> & x, Random * theRandom);

  // to cumulative (for relSigma)
  void toCumulative();

  //
  void normalizeTrue();

  //
  bool byVol = false;  // divide by volume (for physics)
  bool toWrite = true; //

 private:
  Histo(const Histo &);
  Histo(Histo &);

  float getVolume();

  void add(std::pair<float,float> & a, const float & v);
  void div(std::pair<float,float> & a, const std::pair<float,float> & b);

  bool index(const std::vector<float> & x, std::vector<int> & ix);

  // write, read
  void write(std::ofstream & file, std::vector<int>  & ix,
                                   std::vector<float> & x, int j);
  void write(    ogzstream & file, std::vector<int>  & ix,
                                   std::vector<float> & x, int j);

  void read (    igzstream & file, std::vector<int>  & ix, int j);
  void read (std::ifstream & file, std::vector<int>  & ix, int j);

  // data containers (val,sig2)
                                      std::vector<std::pair<float,float>>    a1;
                          std::vector<std::vector<std::pair<float,float>>>   a2;
              std::vector<std::vector<std::vector<std::pair<float,float>>>>  a3;
  std::vector<std::vector<std::vector<std::vector<std::pair<float,float>>>>> a4;
  std::vector<std::vector<std::vector<std::vector<        unsigned char >>>> u4;

  int type;
  std::string name;
  bool norm;
  bool zip;
  float binvolume;
  float sum; // sum of line, if norm

  std::mutex mtx;
};

#endif

