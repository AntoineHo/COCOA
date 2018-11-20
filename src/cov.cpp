#include "cov.h"

// Classes functions

/***
 *       _____ ____  _      ____  _____   _____
 *      / ____/ __ \| |    / __ \|  __ \ / ____|
 *     | |   | |  | | |   | |  | | |__) | (___
 *     | |   | |  | | |   | |  | |  _  / \___ \
 *     | |___| |__| | |___| |__| | | \ \ ____) |
 *      \_____\____/|______\____/|_|  \_\_____/
 *
 *
 */

Color::Color() : r(0.0), g(0.0), b(0.0)
{ }
Color::Color(int _r, int _g, int _b) : r((float) _r / 255), g((float) _g / 255), b((float) _b / 255)
{ }

/***
 *      _____          _____  _       _
 *     |  __ \        |  __ \| |     | |
 *     | |__) | __ ___| |__) | | ___ | |_
 *     |  ___/ '__/ _ \  ___/| |/ _ \| __|
 *     | |   | | |  __/ |    | | (_) | |_
 *     |_|   |_|  \___|_|    |_|\___/ \__|
 *
 *
 */

PrePlot::PrePlot(std::string pname, unsigned int pstart, unsigned int pend) : name(pname), start(pstart), end(pend)
{ }

std::string PrePlot::getName() { return name; }

unsigned int PrePlot::getStart() { return start; }

unsigned int PrePlot::getEnd() { return end; }

/***
 *      _    _ _       _     _      _       _     _
 *     | |  | (_)     | |   | |    (_)     | |   | |
 *     | |__| |_  __ _| |__ | |     _  __ _| |__ | |_
 *     |  __  | |/ _` | '_ \| |    | |/ _` | '_ \| __|
 *     | |  | | | (_| | | | | |____| | (_| | | | | |_
 *     |_|  |_|_|\__, |_| |_|______|_|\__, |_| |_|\__|
 *                __/ |                __/ |
 *               |___/                |___/
 */

Highlight::Highlight() : refName("None"), geneName("None"), start(0), stop(0), col(Color())
{ }
Highlight::Highlight(std::string _refName, std::string _geneName, unsigned int _start, unsigned int _stop, Color _col) : refName(_refName), geneName(_geneName), start(_start), stop(_stop), col(_col)
{ }

/***
 *       _____
 *      / ____|
 *     | |     _____   __
 *     | |    / _ \ \ / /
 *     | |___| (_) \ V /
 *      \_____\___/ \_/
 *
 *
 */

Cov::Cov() : coverage({0}), refName("None"), bamFile("None"), length(0), max(0.0), min(0.0), mean(0.0), stdv(0.0)
{}
Cov::Cov(std::vector<unsigned int> _coverage, std::string _bamFile, std::string _refName, unsigned int _length, double _max, double _min, double _mean, double _stdv) : coverage(_coverage), refName(_refName), bamFile(_bamFile), length(_length), max(_max), min(_min), mean(_mean), stdv(_stdv)
{}
std::vector<unsigned int> Cov::getCov() { return coverage; }
std::string Cov::getBam() { return bamFile; }
std::string Cov::getName() { return refName; }
unsigned int Cov::getLength() { return length; }
double Cov::getMax() { return max; }
double Cov::getMin() { return min; }
double Cov::getMean() { return mean; }
double Cov::getStdv() { return stdv; }

/***
 *       _____  _____   _
 *      / ____|/ ____| | |
 *     | |  __| |      | | __ _ _   _  ___ _ __
 *     | | |_ | |      | |/ _` | | | |/ _ \ '__|
 *     | |__| | |____  | | (_| | |_| |  __/ |
 *      \_____|\_____| |_|\__,_|\__, |\___|_|
 *                               __/ |
 *                              |___/
 */

GClayer::GClayer() : coverage({0}), length(1)
{}
GClayer::GClayer(std::vector<float> _coverage, unsigned int _length) : coverage(_coverage), length(_length)
{}
std::vector<float> GClayer::getGCcov() { return coverage; }
unsigned int GClayer::getLength() { return length; }

/***
 *      _                             _
 *     | |                           | |
 *     | |     __ _ _   _  ___  _   _| |_
 *     | |    / _` | | | |/ _ \| | | | __|
 *     | |___| (_| | |_| | (_) | |_| | |_
 *     |______\__,_|\__, |\___/ \__,_|\__|
 *                   __/ |
 *                  |___/
 */

Layout::Layout() : WIDTH(10200), HEIGHT(2100), refName("None"), Layers(std::vector<Cov> {Cov()}), hLights(std::vector<Highlight> {Highlight()}), gLayer(GClayer())
{ }

Layout::Layout(std::string _refName, std::vector<Cov> _Layers) {
  // VARIABLES FROM ARGS
  refName = _refName;
  Layers = _Layers;
  hLights = std::vector<Highlight> {Highlight()};
  gLayer = GClayer();
  binSize = computeBinsize();
  maxMean = getMaxMean();
  regionLength = Layers[0].getLength();
  palette = {
    Color(204, 0, 0), Color(0, 204, 0), Color(0, 0, 204),
    Color(204,102,0), Color(0,204,102), Color(102,0,204),
    Color(204,204,0), Color(0,204,204), Color(204,0,204),
    Color(102,204,0), Color(0,102,204), Color(204,0,102)
  };
  WIDTH = computeWidth();
  HEIGHT = 2100;
}

Layout::Layout(std::string _refName, std::vector<Cov> _Layers, GClayer _gLayer) {
  // VARIABLES FROM ARGS
  refName = _refName;
  Layers = _Layers;
  hLights = std::vector<Highlight> {Highlight()};
  gLayer = _gLayer;
  binSize = computeBinsize();
  maxMean = getMaxMean();
  regionLength = Layers[0].getLength();
  palette = {
    Color(204, 0, 0), Color(0, 204, 0), Color(0, 0, 204),
    Color(204,102,0), Color(0,204,102), Color(102,0,204),
    Color(204,204,0), Color(0,204,204), Color(204,0,204),
    Color(102,204,0), Color(0,102,204), Color(204,0,102)
  };
  WIDTH = computeWidth();
  HEIGHT = 2100;
}

Layout::Layout(std::string _refName, std::vector<Cov> _Layers, std::vector<Highlight> _hLights) {
  // VARIABLES FROM ARGS
  refName = _refName;
  Layers = _Layers;
  hLights = _hLights;
  gLayer = GClayer();
  binSize = computeBinsize();
  maxMean = getMaxMean();
  regionLength = Layers[0].getLength();
  palette = {
    Color(204, 0, 0), Color(0, 204, 0), Color(0, 0, 204),
    Color(204,102,0), Color(0,204,102), Color(102,0,204),
    Color(204,204,0), Color(0,204,204), Color(204,0,204),
    Color(102,204,0), Color(0,102,204), Color(204,0,102)
  };
  WIDTH = computeWidth();
  HEIGHT = 2100;
}

Layout::Layout(std::string _refName, std::vector<Cov> _Layers, std::vector<Highlight> _hLights, GClayer _gLayer) {
  // VARIABLES FROM ARGS
  refName = _refName;
  Layers = _Layers;
  hLights = _hLights;
  gLayer = _gLayer;
  binSize = computeBinsize();
  maxMean = getMaxMean();
  regionLength = Layers[0].getLength();
  palette = {
    Color(204, 0, 0), Color(0, 204, 0), Color(0, 0, 204),
    Color(204,102,0), Color(0,204,102), Color(102,0,204),
    Color(204,204,0), Color(0,204,204), Color(204,0,204),
    Color(102,204,0), Color(0,102,204), Color(204,0,102)
  };
  WIDTH = computeWidth();
  HEIGHT = 2100;
}

void Layout::plotInit(unsigned int begin, unsigned int end) {

  surface = cairo_svg_surface_create(std::string(refName + "_" + toString(begin) + "-" + toString(end) + ".cov.svg").c_str(), WIDTH, HEIGHT);
  cr = cairo_create(surface);

  cairo_set_line_width(cr, 5);                      // SET LINE WIDTH
  cairo_rectangle(cr, 0, 0, WIDTH, HEIGHT);         // RECTANGLE BACKGROUND
  cairo_set_source_rgb(cr, 1, 1, 1);                // SETS WHITE COLOR
  cairo_fill(cr);                                   // FILLS THE RECTANGLE

  cairo_set_source_rgb(cr, 0, 0, 0);                // RESETS TO BLACK
  cairo_move_to(cr, 95, 1902.5);                    // ORIGIN ALWAYS AT 95 - 1902.5
  cairo_line_to(cr, WIDTH - 97.5, 1902.5);          // X END = WIDTH - 100 + 2.5
  cairo_stroke(cr);                                 // PRINT LINE
  cairo_move_to(cr, 97.5, 1902.5);                  // RESETS BRUSH TO ORIGIN
  cairo_line_to(cr, 97.5, 97.5);                    // Y END ALWAYS AT 97.5
  cairo_stroke(cr);                                 // PRINT LINE

  cairo_rectangle(cr, 100, 1940, WIDTH - 200, 55);  // PLOT REF (LIGHT GREY CHROM)
  //cairo_set_source_rgb(cr, 0.85, 0.85, 0.85);       // SETS LIGHTGREY
  //cairo_fill_preserve(cr);                          // FILL LIGHTGREY
  cairo_set_source_rgb(cr, 0.25, 0.25, 0.25);       // SETS DARKGREY
  cairo_stroke(cr);                                 // STROKE DARKGREY

  std::string yAxLabel("Coverage");                 // SETS Y AX TEXT
  cairo_set_font_size (cr, 72.0);                   // CHANGE FONT SIZE
  cairo_rotate(cr, 90);                             // ROTATES THE TEXT
  cairo_move_to (cr, 1000, 1000);                   // POSITION OF Y AX LABEL
  cairo_rotate(cr, -90);                            // ROTATE BACK
  cairo_show_text (cr, yAxLabel.c_str());           // SHOW TEXT
  cairo_move_to (cr, WIDTH / 2 , 2058);             // POSITION OF X AX LABEL
  cairo_show_text (cr, refName.c_str());            // SHOW TEXT

  cairo_set_font_size (cr, 18.0);                   // X TICKS LABEL FONTSIZE
  for (int xp(100); xp <= WIDTH - 100; xp+=100) {   // FOR EACH 100 BINS SET POSITION
    cairo_set_source_rgb(cr, 0, 0, 0);              // CHANGE COLOR TO BLACK
    cairo_move_to(cr, xp, 1902.5);                  // MOVES TO Y POSITION OF X TICK
    cairo_line_to(cr, xp, 1912.5);                  // PLOT X TICK LINE
    cairo_stroke(cr);                               // STROKE LINE
    std::string xTickLabel(toString(xp-100));       // SET LABEL
    cairo_text_extents(cr, xTickLabel.c_str(), &extents);           // SET TEXT EXTENT
    double xTextCenter = xp-(extents.width/2 + extents.x_bearing);  // CENTER THE X TICK TEXT
    cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);        // SETS COLOR TO DARKGREY
    cairo_move_to (cr, xTextCenter, 1927.5);        // SETS POSITION OF X TICK LABEL
    cairo_show_text (cr, xTickLabel.c_str());       // SHOW LABEL
  }

  unsigned int maxYscale( (unsigned int) 1.8 * maxMean ); // Mean + 100% (so that most of the dots are in the middle part of the y axis graph)
  while (maxYscale % 36 != 0) { maxYscale++; }
  unsigned int tickStep(maxYscale/36);      // SETS Y AX TICKS STEP
  if (maxYscale < 36) { maxYscale = 36; }   // JUST IN CASE MAX IS < 36
  if (tickStep < 1) { tickStep = 1; }       // THEN ALSO RESET TICKSTEP TO INTEGER 1

  cairo_set_font_size (cr, 18.0);           // CHANGE FONT SIZE
  int t(0);                                 // Counts ticks
  for (int yp(100); yp <= 1900; yp+=50) {   // LOOPS 50* POSITIONS BETWEEN 100 & 1900
    cairo_set_source_rgb(cr, 0, 0, 0);      // SETS COLOR TO BLACK
    cairo_move_to(cr, 87.5, yp);            // X TICK POSITION
    cairo_line_to(cr, 97.5, yp);            // X TICK LINE POSITION
    cairo_stroke(cr);                       // STROKE LINE
    std::string yTickLabel(toString(maxYscale-(t*tickStep)));       // SETS LABEL TEXT
    cairo_text_extents(cr, yTickLabel.c_str(), &extents);           // GETS TEXT EXTENT
    double xTextCenter = 80-extents.width;                          // SETS X CENTERING OF TEXT
    double yTextCenter = yp-(extents.height/2 + extents.y_bearing); // SETS Y CENTERING OF TEXT
    cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);                        // CHANGE COLOR
    cairo_move_to (cr, xTextCenter, yTextCenter);                   // X POSITION OF THE TEXT
    cairo_show_text (cr, yTickLabel.c_str());                       // SHOW TICK TEXT
    t++;
  }
}

void Layout::plotGC() {
  // POT THE GC SCALE (ON THE RIGHT)
  cairo_set_line_width(cr, 5);                                      // SETS LINE WIDTH
  cairo_set_source_rgb(cr, 0, 0, 0);                                // SETS COLOR TO BLACK
  cairo_move_to(cr, WIDTH-97.5, 1905);                              // GC AX START
  cairo_line_to(cr, WIDTH-97.5, 897.5);                             // GC AX END (GC AX IS OF LENGTH 1000)
  // TICKS
  cairo_set_font_size (cr, 18.0);                                   // CHANGE FONT SIZE
  int t(0);                                                         // Counts ticks
  for (int yp(900); yp <= 1900; yp+=100) {                          // LOOPS 50* POSITIONS BETWEEN 100 & 1900
    cairo_set_source_rgb(cr, 0, 0, 0);                              // SETS COLOR TO BLACK
    cairo_move_to(cr, WIDTH-95, yp);                                // X TICK POSITION
    cairo_line_to(cr, WIDTH-85, yp);                                // X TICK LINE POSITION
    cairo_stroke(cr);                                               // STROKE LINE
    std::string yTickLabel(toString(1.0-(t*0.1)));                  // SETS LABEL TEXT
    cairo_text_extents(cr, yTickLabel.c_str(), &extents);           // GETS TEXT EXTENT
    double xTextCenter = WIDTH-80;                                  // SETS X CENTERING OF TEXT
    double yTextCenter = yp-(extents.height/2 + extents.y_bearing); // SETS Y CENTERING OF TEXT
    cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);                        // CHANGE COLOR
    cairo_move_to (cr, xTextCenter, yTextCenter);                   // X POSITION OF THE TEXT
    cairo_show_text (cr, yTickLabel.c_str());                       // SHOW TICK TEXT
    t++;
  }

  // GC CONTENT
  cairo_set_source_rgba(cr, 0.3, 0.3, 0.3, 0.75);                   // Solid color
  cairo_set_line_width(cr, 1);

  std::vector<float> gc(gLayer.getGCcov());
  unsigned int gSize( gc.size() );

  if ( binSize == 1.0 ) {
    for (unsigned int p(0); p < gSize; p++) {
      if (p == 0) { cairo_move_to( cr, p+100, compGCYpos(gc[p]) ); }
      cairo_line_to(cr, p+100, compGCYpos(gc[p]));
    }
    cairo_stroke(cr);
  } else if ( binSize != 1.0 ) {
    std::vector<std::vector<float>> bins( splitVector(gc, binSize) ); // binSize is computed at init of Layout object
    std::vector<float> mVec;                                          // INITS A MEAN VECTOR
    for (int j(0); j < bins.size(); j++ ) {                           // FOR EACH BIN
      double mean( meanV(bins[j]) );                                  // COMPUTE MEAN
      mVec.push_back(mean);                                           // PUSH BACK MEAN TO A VECTOR
    }
    unsigned int vecSize(mVec.size());
    for (unsigned int p(0); p < vecSize; p++) {
      if (p == 0) { cairo_move_to( cr, p+100, compGCYpos(gc[p]) ); }
      cairo_line_to(cr, p+100, compGCYpos( gc[p] ));
    }
    cairo_stroke(cr);
  }
}

int Layout::compGCYpos(float gc) {
  return 1900 - (gc*(900));
}

void Layout::plotHglt() {
  // MUST CHECK IF REGION IS IN SPAN OF REGION TO PLOT
  unsigned int hSize(hLights.size());             // GET HIGHLIGHT VECTOR SIZE
  Highlight cHglt;
  unsigned int refSize(Layers[0].getLength());    // Gets size of reference
  if ( binSize == 1.0 ) {
    for (unsigned int p(0); p < hSize; p++) {
      cHglt = hLights[p];
      cairo_set_source_rgba(cr, cHglt.col.r, cHglt.col.g, cHglt.col.b, 0.5);
      float rS(cHglt.start);
      float rLen(cHglt.stop-rS);
      cairo_rectangle(cr, rS + 100, 1942.5, rLen, 50); // PLOT REF (LIGHT GREY CHROM)
      cairo_fill_preserve(cr);                         // FILL LIGHTGREY
      cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);         // SETS DARKGREY
      cairo_stroke(cr);                                // STROKE DARKGREY
    }
  }
  else if ( binSize != 1.0 ) {
    for (unsigned int p(0); p < hSize; p++) {
      cHglt = hLights[p];
      float rS( ((float) cHglt.start/binSize) );
      float rLen( ((float) cHglt.stop/binSize) - rS );
      cairo_set_source_rgba(cr, cHglt.col.r, cHglt.col.g, cHglt.col.b, 0.5);
      cairo_rectangle(cr, rS + 100, 1942.5, rLen, 50);  // PLOT REF (LIGHT GREY CHROM)
      cairo_fill_preserve(cr);                                        // FILL LIGHTGREY
      cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);                        // SETS DARKGREY
      cairo_stroke(cr);                                               // STROKE DARKGREY
    }
  }
}

void Layout::plotLayers() {
  // FOR EACH LAYER IN LAYERS
  Color cCol;                                             // INITS A COLOR OBJECT
  unsigned int maxYscale( (unsigned int) 2 * maxMean );
  for (int i(0); i < Layers.size(); i++) {
    std::vector<unsigned int> vCov(Layers[i].getCov()); // GET COVERAGE INFO
    if ( i < palette.size() ) {
      cCol = palette[i];     // GETS THE COLOR FROM PALETTE
    } else {
      cCol = Color(200,200,200);
    }
    // CREATE bins & mean vector if length > 10000
    unsigned int covSive(vCov.size());
    if ( binSize == 1.0 ) {
      plotBinsize( binSize );
      cairo_set_line_width(cr, 1);
      int cCov(0);
      for (int p(0); p < covSive; p++) {
        cCov = vCov[p];
        if (cCov <= maxYscale) {
          cairo_set_source_rgba(cr, cCol.r, cCol.g, cCol.b, 0.5);  // Solid color
          cairo_move_to(cr, p+100, 1900);
          cairo_line_to(cr, p+100, compYpos(cCov, maxYscale));
          cairo_stroke(cr);
        } else {
          cairo_set_source_rgba(cr, 0.1, 0.1, 0.1, 0.6);  // Solid color
          cairo_move_to(cr, p+100, 1900);
          cairo_line_to(cr, p+100, 100);
          cairo_stroke(cr);
        }
      }
    } else if ( binSize != 1.0) {
      plotBinsize( binSize );
      std::vector<std::vector<unsigned int>> bins( splitVector(vCov, binSize) ); // binSize is computed at init of Layout object
      std::vector<double> mVec;                 // INITS A MEAN VECTOR
      for (int j(0); j < bins.size(); j++ ) {   // FOR EACH BIN
        double mean( meanV(bins[j]) );          // COMPUTE MEAN
        mVec.push_back(mean);                   // PUSH BACK MEAN TO A VECTOR
      }
      unsigned int vecSize(mVec.size());
      cairo_set_line_width(cr, 1);
      int cCov(0);
      for (int p(0); p < vecSize; p++) {
        cCov = mVec[p];
        if (cCov <= maxYscale) {
          cairo_set_source_rgba(cr, cCol.r, cCol.g, cCol.b, 0.5);  // Solid color
          cairo_move_to(cr, p+100, 1900);
          cairo_line_to(cr, p+100, compYpos(cCov, maxYscale));
          cairo_stroke(cr);
        } else {
          cairo_set_source_rgba(cr, 0.1, 0.1, 0.1, 0.6);  // Solid color
          cairo_move_to(cr, p+100, 1900);
          cairo_line_to(cr, p+100, 100);
          cairo_stroke(cr);
        }
      }
    }
  }
}

void Layout::plotOut() {
  cairo_destroy(cr);
  cairo_surface_destroy(surface);
}

void Layout::plotBinsize(float binFloat) {
  cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);                // SETS COLOR TO STRONG DARKGREY
  cairo_select_font_face (cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size (cr, 50.0);                         // SETS THE FONT SIZE
  std::string binLabel("binsize: " + toString(binFloat)); // SETS THE LABEL
  cairo_move_to (cr, 20, 70);                             // SETS POSITION OF LABEL
  cairo_show_text (cr, binLabel.c_str());                 // SHOW TEXT
}

float Layout::computeBinsize() {
  float i(1.0);
  unsigned int length(Layers[0].getLength());
  while ( (length / i) > 10000 ) { i++; } // Here I put 20 000 because minimal bin size is 2 if binsize is lower than 2 then I will just plot more points
  std::cout << "\tAUTO: bin size for coverage = " << i << "bp" << std::endl;
  return i;
}

double Layout::getMaxMean() {
  double maxMean(0);
  for (int i(0); i < Layers.size(); i++) {
    double cMean( Layers[i].getMean() );
    if ( cMean > maxMean ) { maxMean = cMean; }
  }
  return maxMean;
}

std::vector<double> Layout::binVector(Cov& oCov, int binsize) {
  std::vector<unsigned int> vCov( oCov.getCov() );                            // GETS COVERAGE VECTOR
  std::vector<std::vector<unsigned int>> bins( splitVector(vCov, binsize) );  // CREATES BINS
  std::vector<double> mVec;                                                   // INITS FINAL VECTOR OF MEANS
  for (int i(0); i < bins.size(); i++) {                                      // FOR EACH SUB VECTOR
    double cMean( meanV(bins[i]) );                                           // COMPUTE MEAN
    mVec.push_back(cMean);                                                    // Pushes back current coverage mean to the mean vector
  }
  return mVec;                                                                // Returns the vector
}

int Layout::compYpos(unsigned int coverage, unsigned int maxYscale) {
  return (int) (1900-((1800*coverage)/maxYscale));
}

int Layout::computeWidth() {
  unsigned int length(Layers[0].getLength());
  int _WIDTH(((int) length/binSize) + 201);
  return _WIDTH;
}
