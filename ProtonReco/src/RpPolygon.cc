#include "../interface/RpPolygon.h"

#include <iostream>
#include <fstream>
#include <cmath>

#define sqr(x) ((x)*(x))

#undef Debug

using namespace std;

/*****************************************************************************/
RpPolygon::RpPolygon()
{
}

/*****************************************************************************/
RpPolygon::~RpPolygon()
{
}

/*****************************************************************************/
// check if point is below the line
bool RpPolygon::getBelow(const vector<double> & line,
                         const vector<double> & point)
{
  const double & a = line[0];
  const double & b = line[1];

  const double & x = point[0];
  const double & y = point[1];

  return (y < a*x + b);
}

/*****************************************************************************/
// calculate intersection of line segment with straight line
vector<double> RpPolygon::getIntersect(const vector<vector<double>> & segment,
                                       const vector<double> & line)
{
  const double & P1x = segment[0][0];
  const double & P1y = segment[0][1];
  const double & P2x = segment[1][0];
  const double & P2y = segment[1][1];

  const double & a = line[0];
  const double & b = line[1];

  double lam = (a*P1x+b - P1y)/((P2y - P1y) - a*(P2x - P1x));

  const double eps = 1e-6;
  if(lam >= 0-eps  && lam <= 1+eps)
  {
    double x = P1x + lam*(P2x - P1x);
    double y = P1y + lam*(P2y - P1y);

    return {x,y};
  }
  else
  {
    cerr << " problem with getIntersect (lam = " << lam << ")" << endl;
    exit(1);
  }
}

/*****************************************************************************/
// print closed polygon
void RpPolygon::printPolygon(const vector<vector<double>> & polygon)
{
  const int n = polygon.size();

  ofstream file("../out/rp/polygon/poly.dat");

  if(n > 0)
  for(int i = 0; i <= n; i++)
    file << " " << polygon[i%n][0]
         << " " << polygon[i%n][1]
         << " " << 1
         << endl;
  else
  {
    file << " 0 0 0" << endl;
#ifdef Debug
    cerr << " no polygon" << endl;
#endif
  }

  file.close();
}

/*****************************************************************************/
// add a line, find resulting polygon
void RpPolygon::addLine(vector<vector<double>> & polygon,
                        const vector<double> & line, bool beBelow)
{
#ifdef Debug
  cerr << " line " << line[0] << " " << line[1] << " " << beBelow << endl;
#endif

  // check relative locations
  vector<bool> isBelow;
  bool allBelow = true, allAbove = true;

  for(auto & p : polygon)
  {
    bool b = getBelow(line,p);

    if(!b) allBelow = false;
    if( b) allAbove = false;

    isBelow.push_back(b);
  }

  int n = polygon.size();

#ifdef Debug
  cerr << "  allBelow = " << allBelow << endl
       << "  allAbove = " << allAbove << endl;
#endif

  // look at crossed segments
  if(!allBelow && !allAbove)
  for(int i = 0; i < n; i++)
  {
    int i0 =  i    % n;
    int i1 = (i+1) % n;

#ifdef Debug
    cerr << "  looking at segment ("<<i0<<","<<i1<<")" << endl;
#endif

    if(isBelow[i0] ^ isBelow[i1])
    {
#ifdef Debug
      cerr << "  intersect with ("<<i0<<","<<i1<<")" << endl;
#endif

      vector<double> p = getIntersect({polygon[i0],polygon[i1]}, line);

#ifdef Debug
      cerr << "  insert at " << i0+1 << endl;
#endif

      polygon.insert(polygon.begin() + i0+1, p      );
      isBelow.insert(isBelow.begin() + i0+1, beBelow);
      i++;

      n = polygon.size(); // update
    }

#ifdef Debug
    cerr << "  polygon.size() = " << n << endl;
#endif
  }

  // erase some
  for(int i = 0; i < n; i++)
  if(isBelow[i] != beBelow)
  {
#ifdef Debug
    cerr << "  erase " << i << endl;
#endif

    polygon.erase(polygon.begin() + i);
    isBelow.erase(isBelow.begin() + i);

    i--; n--;
  }
}

/*****************************************************************************/
// add band to polygon
void RpPolygon::addBand(vector<vector<double>> & polygon,
                        const vector<double> & band)
{
  const double & a = band[0];
  const double & b = band[1];
  const double & d = band[2];
 
  addLine(polygon, {a,b+d},true ); // below
  addLine(polygon, {a,b-d},false); // above
}

/*****************************************************************************/
// calculate intersection of two straight lines
vector<double> RpPolygon::getIntersect(const vector<double> & line1,
                                       const vector<double> & line2)
{
  const double & a1 = line1[0];
  const double & b1 = line1[1];

  const double & a2 = line2[0];
  const double & b2 = line2[1];

  double x = (b1 - b2)/(a2 - a1);
  double y = a1*x + b1;

  return {x,y};
}

/*****************************************************************************/
// find initial parallelogram using the first two bands
vector<vector<double>> RpPolygon::getParallelogram(
  const vector<vector<double>> & bands)
{
  const vector<double> line1_up = {bands[0][0], bands[0][1] + bands[0][2]};
  const vector<double> line1_do = {bands[0][0], bands[0][1] - bands[0][2]};

  const vector<double> line2_up = {bands[1][0], bands[1][1] + bands[1][2]};
  const vector<double> line2_do = {bands[1][0], bands[1][1] - bands[1][2]};

  vector<vector<double>> polygon;

  polygon.push_back(getIntersect(line1_do,line2_up));
  polygon.push_back(getIntersect(line1_do,line2_do));
  polygon.push_back(getIntersect(line1_up,line2_do));
  polygon.push_back(getIntersect(line1_up,line2_up));

  return polygon;
}

/*****************************************************************************/
// calculate centroid and moment of intertia wrt centroid in y direction
LocalFit RpPolygon::calculate(const vector<vector<double>> & polygon)
{
  const int n = polygon.size();

  double Cx=0, Cy=0, Ix=0, Iy=0;

  if(n > 0)
  {
  // https://en.wikipedia.org/wiki/Centroid#Of_a_polygon

  // area
  double A=0;
  for(int i = 0; i < n; i++)
  {
    int i0 =  i    % n;
    int i1 = (i+1) % n;

    const double & x0 = polygon[i0][0];
    const double & x1 = polygon[i1][0];

    const double & y0 = polygon[i0][1];
    const double & y1 = polygon[i1][1];

    A += 1./2 * (x0*y1 - x1*y0);

    Cx += (x0 + x1)*(x0*y1 - x1*y0);
    Cy += (y0 + y1)*(x0*y1 - x1*y0);
  }

  // normalize
  if(A != 0)
  {
    Cx /= 6*A; Cy /= 6*A;
  } else {
    Cx = polygon[0][0];
    Cy = polygon[0][1];
  }

  // https://en.wikipedia.org/wiki/Second_moment_of_area#Any_polygon
  // moments of intertia
 
  for(int i = 0; i < n; i++)
  {
    int i0 =  i    % n;
    int i1 = (i+1) % n;

    {
      const double & x0 = polygon[i0][0] - Cx; // wrt centroid
      const double & x1 = polygon[i1][0] - Cx; // wrt centroid

      const double & y0 = polygon[i0][1];
      const double & y1 = polygon[i1][1];

      Ix += (y0*x1 - y1*x0) * (sqr(x0) + x0*x1 + sqr(x1));
    }

    {
      const double & x0 = polygon[i0][0];
      const double & x1 = polygon[i1][0];

      const double & y0 = polygon[i0][1] - Cy; // wrt centroid
      const double & y1 = polygon[i1][1] - Cy; // wrt centroid

      Iy += (x0*y1 - x1*y0) * (sqr(y0) + y0*y1 + sqr(y1));
    }
  }

  // normalize
  if(A != 0)
  {
    Ix /= 12*A; Iy /= 12*A;
  } else {
    Ix = 0; Iy = 0; // Vx = Vy = 0
  }
  }
//  else cerr << " no polygon" << endl;
//  else n=0, no polygon, no common area of the bands, need for simplex

  //
  return {Cx, Cy, fabs(Ix), fabs(Iy), 0}; // set value to zero
}

/*****************************************************************************/
// find common area, return Cx,Cy,Vx,Vy
LocalFit RpPolygon::findPolygon(const vector<vector<double>> & bands,
                                      bool print)
{
  vector<vector<double>> polygon = getParallelogram(bands);

  for(size_t i = 2; i < bands.size(); i++)
    addBand(polygon, bands[i]);

  // Cx, Cy, sigma_y^2
  LocalFit C = calculate(polygon);

  if(print)
    printPolygon(polygon);

  return C;
}

