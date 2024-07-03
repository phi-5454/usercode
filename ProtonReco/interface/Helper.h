#ifndef _Helper_h_
#define _Helper_h_

#include <vector>
#include <string>

//
class Helper
{
 public:
  Helper();
  ~Helper();

  static std::vector<std::string> parse(std::string line,
                                  const std::string & delimiter);

  static void propeller(int i, int bin);

  static std::string col(int c);
  static std::string col();

  static std::vector<std::vector<std::string>>
    getStreams(int nCpus, bool paraOnly);

 private:
};

#endif
