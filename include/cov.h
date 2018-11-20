#ifndef COV_H_INCLUDED
#define COV_H_INCLUDED

// STANDARD MODULES
#include <iostream>   // Necessary to read/write from/to stdin/out
#include <string>
#include <vector>
#include <cmath>
#include <sstream>    // Necessary to use std::to_string
#include <iterator> // Necessary to read & fill from cin // std::begin, std::end
#include <algorithm> // std::min_element
#include <numeric> // Necessary to use accumulate

// CAIRO Library
#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>

// OBJECTS
class Color {
public:
  Color();
  Color(int r, int g, int b);
  float r;
  float g;
  float b;
};

class PrePlot {
public:
  PrePlot(std::string name, unsigned int start, unsigned int end);
  std::string getName(); // Name of the reference to plot
  unsigned int getStart();
  unsigned int getEnd();
private:
  std::string name;
  unsigned int start;
  unsigned int end;
};

class Highlight {
public:
  Highlight();
  Highlight(std::string refName, std::string geneName, unsigned int start, unsigned int stop, Color col = Color());
  std::string geneName;
  std::string refName;
  unsigned int start;     // ON THE REGION DEFINED BY BED FILE!!!
  unsigned int stop;      // ON THE REGION
  Color col;
};

class Cov {
public:
  Cov();
  Cov(std::vector<unsigned int> coverage, std::string bamFile, std::string refName, unsigned int length, double max, double min, double mean, double stdv);
  std::vector<unsigned int> getCov();
  std::string getBam();
  std::string getName();
  unsigned int getLength();
  double getMax();
  double getMin();
  double getMean();
  double getStdv();
private:
  std::vector<unsigned int> coverage;
  std::string bamFile;
  std::string refName;
  unsigned int length;
  double max;
  double min;
  double mean;
  double stdv;
};

class GClayer {
public:
  GClayer();
  GClayer(std::vector<float> coverage, unsigned int length);
  std::vector<float> getGCcov();
  unsigned int getLength();
private:
  std::vector<float> coverage;
  unsigned int length;
};

class Layout {
public:
  Layout();
  Layout(std::string refName, std::vector<Cov> Layers);
  Layout(std::string refName, std::vector<Cov> Layers, GClayer gLayer);
  Layout(std::string refName, std::vector<Cov> Layers, std::vector<Highlight> hLights);
  Layout(std::string refName, std::vector<Cov> Layers, std::vector<Highlight> hLights, GClayer gLayer);
  std::vector<double> binVector(Cov& oCov, int binsize);
  float computeBinsize();
  double getMaxMean();
  void plotInit(unsigned int begin, unsigned int end);
  void plotGC();
  void plotHglt();
  void plotLayers();
  void plotOut();
  void plotBinsize(float binFloat);
  int compYpos(unsigned int coverage, unsigned int maxYscale);
  int compGCYpos(float gc);
  int computeWidth();
private:
  cairo_surface_t *surface;
  cairo_t *cr;
  cairo_text_extents_t extents;
  int WIDTH;
  int HEIGHT;
  std::string refName;
  std::vector<Cov> Layers;
  std::vector<Highlight> hLights;
  GClayer gLayer;
  std::vector<Color> palette;
  float binSize;
  double maxMean;
  unsigned int regionLength;
};

/***
 *      _______                   _       _
 *     |__   __|                 | |     | |
 *        | | ___ _ __ ___  _ __ | | __ _| |_ ___  ___
 *        | |/ _ \ '_ ` _ \| '_ \| |/ _` | __/ _ \/ __|
 *        | |  __/ | | | | | |_) | | (_| | ||  __/\__ \
 *        |_|\___|_| |_| |_| .__/|_|\__,_|\__\___||___/
 *                         | |
 *                         |_|
 */

template<typename T>
std::vector<std::vector<T>> splitVector(const std::vector<T>& vec, size_t n) {
  // Inits a vector of vectors of type
  std::vector<std::vector<T>> outVec;
  size_t size( (vec.size() - 1) / n + 1 );
  size_t remain( vec.size() % n );

  for ( size_t i(0); i < size; i++ ) {
    std::vector<T> bin;
    auto start = std::next(vec.begin(), i*n);
    auto end = std::next(vec.cbegin(), i*n + n);
    std::copy( start, end, std::back_inserter(bin) );
    outVec.push_back(bin);
  }

  return outVec;
}

template<typename T>
T meanV(const std::vector<T>& vec) {
  T sum(0);
  int size(vec.size());
  for(int i(0); i < size; ++i) {
    sum += vec[i];
  }
  return sum/size;
}

template<typename T>
std::string toString(T toLine) {
  std::string toReturn;
  std::stringstream strstream;
  strstream << toLine;
  strstream >> toReturn;
  return toReturn;
}

#endif
