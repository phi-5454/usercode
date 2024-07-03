#include "../interface/Helper.h"

#include <iostream>
#include <sstream>

using namespace std;

/*****************************************************************************/
Helper::Helper()
{
}

/*****************************************************************************/
Helper::~Helper()
{
}

/*****************************************************************************/
vector<string> Helper::parse(string line, const string & delimiter)
{
  vector<string> list;

  size_t pos = 0;
  while((pos = line.find(delimiter)) != string::npos)
  {
    string token = line.substr(0, pos);

    if(token != "")
      list.push_back(token);

    line.erase(0, pos + delimiter.length());
  }

  list.push_back(line); // remainder

  return list;
}

/*****************************************************************************/
void Helper::propeller(int i, int big)
{
  char c[4] = {'|', '/','-','\\'};

  int small = big/20;

  if(i % small == 0) cerr << c[(i / small) % 4] << "\b";
  if(i % big   == 0) cerr << ".";
}

/*****************************************************************************/
string Helper::col(int c)
{
  ostringstream os;
  os << "\033[3" << c << "m";
  return os.str();
}

string Helper::col()
{
  ostringstream os;
  os << "\033[0m";
  return os.str();
}

/*****************************************************************************/
vector<vector<string>> Helper::getStreams(int nCpus, bool paraOnly)
{
 // setup lists
  vector<vector<string>> lists;

  if(paraOnly)
  { // process TT and BB, for veto calculation
    if(nCpus >= 4)
      lists = { {"TOTEM40"},{"TOTEM41"},{"TOTEM42"},{"TOTEM43"} };

    if(nCpus == 2)
      lists = { {"TOTEM40",  "TOTEM41"},{"TOTEM42",  "TOTEM43"} };

    if(nCpus == 1)
      lists = { {"TOTEM40",  "TOTEM41",  "TOTEM42",  "TOTEM43"} };
  }
  else
  { // process all
    if(nCpus >= 4)
      lists = { {"TOTEM20","TOTEM40"}, {"TOTEM21","TOTEM41"},
                {"TOTEM22","TOTEM42"}, {"TOTEM23","TOTEM43"} };

    if(nCpus == 2)
      lists = { {"TOTEM20","TOTEM40",   "TOTEM21","TOTEM41"},
                {"TOTEM22","TOTEM42",   "TOTEM23","TOTEM43"} };

    if(nCpus == 1)
      lists = { {"TOTEM20","TOTEM40",   "TOTEM21","TOTEM41",
                 "TOTEM22","TOTEM42",   "TOTEM23","TOTEM43"} };

    if(nCpus == 0) // 2print
      lists = { {"TOTEM20","TOTEM40"} };
  }

  return lists;
}

