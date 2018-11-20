#include "opts.h"

// Classes functions
OOpts::OOpts(std::vector<std::string> _bams, std::string _bed, std::string _highlight, std::string _reference, int _phred, bool _GC, bool _Q, bool _H, bool _write, bool _plot) : bams(_bams), bed(_bed), highlight(_highlight), reference(_reference), phred(_phred), GC(_GC), Q(_Q), H(_H), write(_write), plot(_plot)
{
  if (reference != "") { hasref = true; } else { hasref = false; }
}

std::vector<std::string> OOpts::getBAMs() { return bams; }
std::string OOpts::getBED() { return bed; }
std::string OOpts::getHighlight() { return highlight; }
std::string OOpts::getRef() { return reference; }
int OOpts::getMinPhred() { return phred; }
bool OOpts::plotGC() { return GC; }
bool OOpts::plotQual() { return Q; }
bool OOpts::plotHighlight() { return H; }
bool OOpts::doPlot() { return plot; }
bool OOpts::hasRef() { return hasref; }
bool OOpts::writeToFile() { return write; }

// FUNCTIONS
OOpts readOpts(int argc, char* argv[]) {
  // REQUIRED BOOLS
  bool gotBAM(false);
  bool gotBED(false);
  std::vector<std::string> bams;
  std::string bed("");
  std::string highlight("");
  std::string reference("");
  bool plotGC(false);
  bool plotQual(false);
  bool plotHighlight(false);
  bool gotReference(false);
  bool writeToFile(false);
  int phred(13);
  for (int i(0); i < argc; i++) {
    // GET BAMS
    if (std::string(argv[i]) == "--bam" || std::string(argv[i]) == "-bm") {
      // COMMA SEPARATED BAM FILES
      i++;
      int p(0);
      char *tkn(strtok(argv[i], ","));
      while (tkn != NULL) {
        bams.push_back(std::string(tkn));
        p++;
        tkn = strtok(NULL, ",");
      }
      if ( bams.size() > 0 ) { gotBAM = true; } // Checks if bams was filled
    } else if (std::string(argv[i]) == "--bed" || std::string(argv[i]) == "-bd") {
      i++;
      bed = std::string(argv[i]);
      if ( bed != "" ) { gotBED = true; } // Checks if bed was found
    } else if (std::string(argv[i]) == "--highlight" || std::string(argv[i]) == "-hl") {
      i++;
      highlight = std::string(argv[i]);
      if ( highlight != "" ) {
        plotHighlight = true;
        OPTS::HIGHLIGHT = true;
      } // Checks if highlight was found
    }else if (std::string(argv[i]) == "--reference" || std::string(argv[i]) == "-r") {
      i++;
      reference = std::string(argv[i]);
      if ( reference != "" ) {
        gotReference = true;
        OPTS::HASREF = true;
      } // Checks if highlight was found
    } else if (std::string(argv[i]) == "--GC" || std::string(argv[i]) == "-g") {
      plotGC = true;
      OPTS::GC = true;
    } else if (std::string(argv[i]) == "--Quality" || std::string(argv[i]) == "-q") {
      plotQual = true;
      OPTS::QUAL = true;
    } else if (std::string(argv[i]) == "--write" || std::string(argv[i]) == "-w") {
      writeToFile = true;
      OPTS::WRITE = true;
    } else if (std::string(argv[i]) == "--phred" || std::string(argv[i]) == "-p") {
      i++;
      phred = strtol(argv[i],NULL,10);
      OPTS::PHRED = phred;
    }
  }
  if ( !gotBAM ) {
    std::cout << "WARNING: No bam file was set!" << std::endl;
    return OOpts(bams, bed, highlight, reference, phred, plotGC, plotQual, plotHighlight, writeToFile, false);
  } else if ( !gotBED ) {
    std::cout << "WARNING: No bed file was set!" << std::endl;
    return OOpts(bams, bed, highlight, reference, phred, plotGC, plotQual, plotHighlight, writeToFile, false);
  } else if ( plotGC && !gotReference ) {
    std::cout << "WARNING: Asked to plot GC without reference genome!" << std::endl;
    return OOpts(bams, bed, highlight, reference, phred, plotGC, plotQual, plotHighlight, writeToFile, false);
  } else {
    return OOpts(bams, bed, highlight, reference, phred, plotGC, plotQual, plotHighlight, writeToFile, true);
  }
}
