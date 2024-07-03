#ifndef _SimplexMinimizer_h_
#define _SimplexMinimizer_h_

#include <functional> 
#include <vector>

const double ftol = 1e-6; // PARAMETER

typedef std::vector<std::vector<double>> Simplex;

class SimplexMinimizer
{
 public:
  SimplexMinimizer(int ndim);
  ~SimplexMinimizer();

  std::pair<std::vector<double>,double> minimize(std::vector<double> point,
                               std::vector<double> dels,
                               std::function<double(std::vector<double>)> & func);

 private:
  void get_psum(const Simplex & p, std::vector<double> & psum);

  double amotry(Simplex & p,
                std::vector<double> & y, std::vector<double> & psum,
                const int ihi, const double fac,
                std::function<double(std::vector<double>)> & func);

  std::pair<std::vector<double>,double> minimize(Simplex p,
                std::function<double(std::vector<double>)> & func);

  std::vector<double> y;

  const int ndim;
  int mpts;
};

#endif
