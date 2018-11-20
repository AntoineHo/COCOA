#ifndef BED_H_INCLUDED
#define BED_H_INCLUDED

// STANDARD MODULES
#include <iostream>   // Necessary to read/write from/to stdin/out
#include <string>     // Necessary to use strings
#include <vector>     // Necessary to use dynamically allocated tables
#include <fstream>    // Necessary to read/write from/to files
#include <cstring>    // Necessary to use c_string functions ex: strtok
#include <cstdlib>    // Necessary to use strtol

// SEQAN MODULES
#include <seqan/basic.h>
#include <seqan/bed_io.h>

// LOCAL MODULES
#include "cov.h" // NECESSARY FOR PREPLOT OBJECTS

// PROTOTYPES
/**
 * \brief Function that returns a vector of PrePlot object and fills a map
 * \param const std::string& bedFilePath a reference to the string containing the .bed filepath
 * \param std::map<std::string, int> * pRefLen a pointer to the map of reference name : reference length
 * \return std::vector<PrePlot>
 **/
std::vector<PrePlot> readBed(const std::string& bedFilePath, std::map<std::string, unsigned int> * pRefLen);

#endif
