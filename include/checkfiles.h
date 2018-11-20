#ifndef CHECKFI_H_INCLUDED
#define CHECKFI_H_INCLUDED


// STANDARD MODULES
#include <iostream>   // Necessary to read/write from/to stdin/out
#include <string>     // Necessary to use strings
#include <vector>     // Necessary to use dynamically allocated tables
#include <fstream>    // Necessary to read/write from/to files
#include <cstring>    // Necessary to use c_string functions ex: strtok

// SEQAN MODULES
#include <seqan/basic.h>
#include <seqan/bam_io.h>
#include <seqan/bed_io.h>
#include <seqan/seq_io.h>

// PROTOTYPES OF FUNCTIONS

/**
 * \brief function checking bam files
 * \param std::string& bamF string containing bam file path
 * \return true: OK false: ERRORS
**/
bool checkBAM(std::string& bamF);

/**
 * \brief function checking sam files
 * \param std::string& bams string containing sam file path
 * \return true: OK false: ERRORS
**/
bool checkSAM(std::string& samF);

/**
 * \brief function that checks reference file
 * \param std::string& string containing the reference FQ/FA file path
 * \return true: OK false: ERRORS
**/
bool checkREF(std::string& ref);

/**
 * \brief function that checks bed file
 * \param std::string& string containing the reference bed file path
 * \return true: OK false: ERRORS
**/
bool checkBED(std::string& bed);

/**
 * \brief Function that returns a true if filepath exists
 * \param filepath const std::string& a reference to the string containing the .info filepath
 * \return bool
 **/
inline bool fileExist(const std::string& filepath);
#endif
