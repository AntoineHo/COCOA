#ifndef OPTS_H_INCLUDED
#define OPTS_H_INCLUDED

// STANDARD MODULES
#include <iostream>   // Necessary to read/write from/to stdin/out
#include <string>     // Necessary to use strings
#include <vector>     // Necessary to use dynamically allocated tables
#include <cstring>    // Necessary to use c_string functions ex: strtok

namespace OPTS {
    inline bool GC(false);
    inline bool WRITE(false);
    inline bool HIGHLIGHT(false);
    inline bool HASREF(false);
    inline bool QUAL(false);
    inline int PHRED(13);
}

// OBJECTS
class OOpts {
public:
  OOpts(std::vector<std::string> bams, std::string bed, std::string highlight, std::string reference, int phred, bool GC, bool Q, bool H,  bool write, bool plot);
  std::vector<std::string> getBAMs();
  std::string getBED();
  std::string getHighlight();
  std::string getRef();
  int getMinPhred();
  bool plotGC();
  bool plotQual();
  bool plotHighlight();
  bool doPlot();
  bool hasRef();
  bool writeToFile();
private:
  std::vector<std::string> bams;
  std::string bed;
  std::string highlight;
  std::string reference;
  int phred;
  bool GC;
  bool Q;
  bool H;
  bool plot;
  bool hasref;
  bool write;
};

// PROTOTYPES OF FUNCTIONS

/**
 * \brief function that reads arguments from command line
 * \param int argc count of the arguments
 * \param char* argv[] char array containing the arguments
 * \return OOpts object
**/
OOpts readOpts(int argc, char* argv[]);
#endif
