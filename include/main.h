#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

// STANDARD MODULES
#include <iostream>   // Necessary to read/write from/to stdin/out
#include <string>     // Necessary to use strings
#include <vector>     // Necessary to use dynamically allocated tables
#include <fstream>    // Necessary to read/write from/to files
#include <cstring>    // Necessary to use c_string functions ex: strtok
#include <cstdlib>    // Necessary to use strtol
#include <unistd.h>   // Necessary to get filesystem
#include <stdlib.h>
#include <algorithm>  // std::sort
#include <cmath>      // Necessary for math functions
#include <sstream>    // Necessary to use std::to_string
#include <map>        // Necessary to use dictionaries
#include <exception>


// CAIRO Library
#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>

// SEQAN MODULES
#include <seqan/basic.h>
#include <seqan/bam_io.h>
#include <seqan/bed_io.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

// LOCAL MODULES
#include "bed.h"
#include "cov.h"
#include "opts.h"
#include "checkfiles.h"

// PROTOTYPES OF FUNCTIONS

/**
 * \brief Function that returns a true if filepath exists
 * \param filepath const std::string& a reference to the string containing the .info filepath
 * \return bool
 **/
inline bool fileExist(const std::string& filepath);

/**
 * \brief Function that just prints help to std::cout
 * \return Nothing
 **/
void printHelp();

/**
 * \brief function that returns a vector of string splitted at one regex
 * \param const string& input a reference to a const string to split
 * \param const string& regex a reference to a const string regular expression
 * \return std::vector<std::string>
**/
std::vector<std::string> split(std::string& input, const char * delimiters);

/**
 * \brief function that searches reference, reads sam and plots coverage from a PrePlot object
 * \param std::vector<unsigned int>& vCov a reference to a vector of coverage
 * \param unsigned int& bS a reference to the starting base pair
 * \param (const char * bamFilePath a reference to the bam file path
 * \param std::string& refName a reference to the name of ref to plot from bed file
 * \return false:FAIL or true:PASS
**/
bool writeCov(std::vector<unsigned int>& vCov, unsigned int& bS, std::string& output, std::string& refName);

/**
 * \brief function that reads the header of a bam file and fills the map
 * \param std::string& bamF a reference to a string path of the bam file
 * \param std::map<std::string, unsigned int> * pRefLen a pointer to a map of key=string name of the reference and value=unsigned int length of the reference
 * \param std::vector<PrePlot>& toPlot reference to the vector of string reference names to plot <-- must be found in bam
 * \return false:FAIL or true:PASS
**/
bool readHeader(std::string& bamF, std::map<std::string, unsigned int> * pRefLen, std::vector<std::string>& refNamesToPlot);

/**
 * \brief function that reads an indexed bam file and pushes back a coverage object to the vector of layers
 * \brief this function also outputs a file /path/to/file/bamfilename_referencename.cov
 * \brief the file output has the samtools depth TSV format:   reference    position    coverage
 * \param std::string& bamF a reference to a string path of the bam file
 * \param PrePlot& pPlot a reference to the preplot object to compute coverage for
 * \param unsigned int& refLength a reference to a constant integer reference length
 * \param std::vector<Cov>& layers a reference to the vector of Coverage objects to fill
 * \return false:FAIL or true:PASS
**/
bool computeBAMCoverage(std::string& bamF, PrePlot& pPlot, unsigned int& refLength, std::vector<Cov>& layers);

/**
 * \brief function that reads a sam file and pushes back a coverage object to the vector of layers
 * \brief this function also outputs a file /path/to/file/bamfilename_referencename.cov
 * \brief the file output has the samtools depth TSV format:   reference    position    coverage
 * \param std::string& samF a reference to a string path of the sam file
 * \param PrePlot& pPlot a reference to the preplot object to compute coverage for
 * \param unsigned int& refLength a reference to a constant integer reference length
 * \param std::vector<Cov>& layers a reference to the vector of Coverage objects to fill
 * \return false:FAIL or true:PASS
**/
bool computeSAMCoverage(std::string& bamF, PrePlot& pPlot, unsigned int& refLength, std::vector<Cov>& layers);

/**
 * \brief function that reads a coverage vector and computes min, max, mean and stdv
 * \param const std::vector<unsigned int>& vCov reference to a const object vCov
 * \param double& stats[] reference to an array of double
**/
void computeStats( const std::vector<unsigned int>& vCov, double (&stats)[4] );

/**
 * \brief function that reads a fa/fq file and returns the region considered in preplot object
 * \param std::string& refFile  reference fasta file (origin of the mapping)
 * \param std::string& refName  reference name of the contig (found in bed file)
 * \param unsigned int& pS      start position of region to plot
 * \param unsigned int& pE      end position of region to plot
 * \param std::string& region   Reference to the string region output of this function
 * \return false:FAIL or true:PASS
**/
bool extractSeq(std::string& refFile, std::string& refName, unsigned int& pS, unsigned int& pE, std::string& region);

/**
 * \brief function that computes GC content in a dna string
 * \param std::string& dna  dna string
 * \return a double ratio of GC content
**/
float computeGC(std::string& dna);

/**
 * \brief function that reads the highlight
 * \param std::vector<Highlight>& vH    vector of highlight objects to fill
 * \param std::string& bedF             string containing bed file path
 * \param unsigned int& pS              int start of to plot region
 * \param unsigned int& pE              int end of to plot region
 * \return false:FAIL or true:PASS
**/
bool readHighlight(std::vector<Highlight>& vH, std::string& bedF, unsigned int& pS, unsigned int& pE);
#endif
