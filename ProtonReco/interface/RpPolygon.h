#ifndef _RpPolygon_h_
#define _RpPolygon_h_

#include <vector>
#include <iostream>
#include <cmath>

//
struct LocalFit {
  double Cx;
  double Cy;
  double Vx;
  double Vy;
  double Val;

  void print()
  {
    std::cerr << " Cx = " << Cx << " +- " << sqrt(Vx) << ","
              << " Cy = " << Cy << " +- " << sqrt(Vy) << ","
              << " χ2 = " << Val << std::endl;
  }
};

//
class RpPolygon
{
 public:
  RpPolygon();
  ~RpPolygon();

  LocalFit findPolygon(
    const std::vector<std::vector<double>> & bands, bool print);

 private:
  bool getBelow(const std::vector<double> & line,
                const std::vector<double> & point);

  std::vector<double> getIntersect(
      const std::vector<std::vector<double>> & segment,
                  const std::vector<double> & line);

  void printPolygon(const std::vector<std::vector<double>> & polygon);

  void addLine(std::vector<std::vector<double>> & polygon,
                     const std::vector<double> & line, bool beBelow);

  void addBand(std::vector<std::vector<double>> & polygon,
                     const std::vector<double> & band);

  std::vector<double> getIntersect(const std::vector<double> & line1,
                                   const std::vector<double> & line2);

  std::vector<std::vector<double>> getParallelogram(
    const std::vector<std::vector<double>> & bands);

  LocalFit calculate(const std::vector<std::vector<double>> & polygon);
};

#endif
