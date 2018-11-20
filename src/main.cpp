#include "main.h"

// Arg1: bam file
// Arg2: bed file or contig_name:start_position-end_position
int main(int argc, char* argv[]) {
  if (argv[1] == "-h" || argv[1] == "--help" || argc < 2) {
    printHelp();
    return 1;
  }

  // ARGUMENT PARSE & CHECK
  OOpts opts( readOpts(argc, argv) );
  if ( !opts.doPlot() ) {
    printHelp();
    return 1;
  }

  bool plotGC( opts.plotGC() );
  bool plotQual( opts.plotQual() );

  // FILE CHECK
  std::cout << "# Checking files..." << std::endl;
  std::vector<std::string> bams( opts.getBAMs() );
  for (int i(0); i < bams.size(); i++) {
    if (bams[i].substr(bams[i].find_last_of(".") + 1) == "bam") {
      if ( !checkBAM(bams[i]) ) {
        return 1;
      }
    } else if (bams[i].substr(bams[i].find_last_of(".") + 1) == "sam") {
      if ( !checkSAM(bams[i]) ) {
        return 1;
      } else {
        std::cout << "\tWARNING: You are trying to read a .SAM file! Finding alignments can take much more time than with indexed .BAM files!" << std::endl;
      }
    } else {
      std::cout << "ERROR: Wrong extension! File should be SAM or BAM" << std::endl;
      return 1;
    }

  }

  std::string bed( opts.getBED() );
  if ( !checkBED(bed) ) {
    return 1;
  }

  std::string ref;
  if ( opts.hasRef() ) {
    ref = opts.getRef();
    if ( !checkREF(ref) ) {
      return 1;
    }
  }

  std::string hglt;
  if ( opts.plotHighlight() ) {
    hglt = opts.getHighlight();
    if ( !checkBED(hglt) ) {
      return 1;
    }
  }
  std::cout << "DONE!" << std::endl;

  // PROCESSING PREPLOT OBJECTS
  std::cout << "# Preprocessing regions from bed file..." << std::endl;
  std::map<std::string, unsigned int> * pRefLen(0); // Pointer to a map REFNAME LENGTH
  pRefLen = new std::map<std::string, unsigned int>; // Allocates memory to pRefLen pointer
  // Creates a vector toPlot of preplot objects
  std::vector<PrePlot> toPlot( readBed((const std::string) bed, pRefLen) );
  // Creates a vector of string of reference names to plot
  std::vector<std::string> refNamesToPlot;
  for (int i(0); i < toPlot.size(); i++) { refNamesToPlot.push_back( toPlot[i].getName() ); }
  std::cout << "DONE: Stored pre - plotting information!" << std::endl;

  // GETTING LENGTH INFORMATION
  std::cout << "# Reading headers of bam files..." << std::endl;
  for (int i(0); i < bams.size(); i++) {
    if ( !readHeader(bams[i], pRefLen, refNamesToPlot) ) { return 1; }
  }
  std::cout << "DONE: found " << toPlot.size() << " reference(s) to plot!" << std::endl;

  // FOR EACH PRE PLOT OBJECT
  for (int i(0); i < toPlot.size(); i++) {
    std::cout << "# Processing: " << toPlot[i].getName() << ":" << toPlot[i].getStart() << "-" << toPlot[i].getEnd() << " ..." << std::endl;
    // Sets a coverage object vector called layers --> ORDER OF BAM FILE IS ORDER OF PLOT FIRST IN FIRST OUT !
    std::vector<Cov> layers;
    // Get necessary info
    std::string ctg(toPlot[i].getName());
    unsigned int pS(toPlot[i].getStart());
    unsigned int pE(toPlot[i].getEnd());
    unsigned int refLength( (*pRefLen)[ ctg ] ); // Get length from reference name of preplot object
    // COMPUTE COVERAGE OF REGION IN EVERY BAM FILE --> STORE TO COVERAGE OBJECT --> STORE THIS OBJECT TO LAYERS
    for (int j(0); j < bams.size(); j++) {
      // FOR EACH BAM FILE IN bams list
      if (bams[j].substr(bams[j].find_last_of(".") + 1) == "bam") {
        if ( !computeBAMCoverage(bams[j], toPlot[i], refLength, layers) ) {
          return 1;
        }
      } else if (bams[j].substr(bams[j].find_last_of(".") + 1) == "sam") {
        if ( !computeSAMCoverage(bams[j], toPlot[i], refLength, layers) ) {
          return 1;
        }
      }
    }
    std::cout << "DONE: Coverage layers: " << layers.size() << std::endl;

    bool GCdone(false);
    // GET GC FROM REFERENCE -> TO LAYER
    GClayer *gcLayer(0);
    gcLayer = new GClayer; // Allocates memory to pRefLen pointer
    if ( OPTS::GC && OPTS::HASREF ) {
      std::cout << "# Extracting reference region and GC computing..." << std::endl;
      std::string region; // REGION TO COMPUTE GC
      if ( !extractSeq(ref, ctg, pS, pE, region) ) {
        return 1;
      }
      unsigned int length( pE-pS );
      int windowSize( 0.01*length ); // Length = 10% of total size
      if (windowSize == 0) {
        windowSize = 1;
      }
      std::cout << "\tAUTO: window size for GC content = " << windowSize << "bp" << std::endl;
      std::vector<float> GC;
      unsigned int x(0);
      while (x+windowSize < length) {
        std::string subRegion(region.substr(x, windowSize));
        GC.push_back( computeGC( subRegion ) ); // Computes GC for substring of reference region
        x++;
      }
      GClayer cGC(GC, (unsigned int) GC.size());
      (*gcLayer) = cGC;
      GCdone = true;
      std::cout << "DONE!" << std::endl;
    } else {
      delete gcLayer;
      gcLayer = 0;
    }

    // FETCH THE HIGHLIGHT LAYER
    bool HGLTdone(false);
    std::vector<Highlight> *hLayer(0);
    hLayer = new std::vector<Highlight>;
    if ( OPTS::HIGHLIGHT )  {
      std::cout << "# Reading highlight regions..." << std::endl;
      if (!readHighlight((*hLayer), hglt, pS, pE)) {
        return 1;
      } else {
        HGLTdone = true;
        std::cout << "Done: found " << (*hLayer).size() << " regions to highlight!" << std::endl;
      }
    } else {
      delete hLayer;
      hLayer = 0;
    }

    // Create a LAYOUT IMAGE
    std::cout << "# Building layout..." << std::endl;
    unsigned int length(toPlot[i].getEnd()-toPlot[i].getStart());

    if ( GCdone && HGLTdone ) {
      Layout lyt( toPlot[i].getName(), layers, (*hLayer), (*gcLayer) );
      lyt.plotInit( toPlot[i].getStart(), toPlot[i].getEnd() );
      lyt.plotLayers();
      lyt.plotHglt();
      lyt.plotGC();
      lyt.plotOut();
    } else if ( GCdone && !HGLTdone ) {
      Layout lyt( toPlot[i].getName(), layers, (*gcLayer) );
      lyt.plotInit( toPlot[i].getStart(), toPlot[i].getEnd() );
      lyt.plotLayers();
      lyt.plotGC();
      lyt.plotOut();
    } else if ( !GCdone && HGLTdone ) {
      Layout lyt( toPlot[i].getName(), layers, (*hLayer) );
      lyt.plotInit( toPlot[i].getStart(), toPlot[i].getEnd() );
      lyt.plotLayers();
      lyt.plotHglt();
      lyt.plotOut();
    } else {
      Layout lyt( toPlot[i].getName(), layers );
      lyt.plotInit( toPlot[i].getStart(), toPlot[i].getEnd() );
      lyt.plotLayers();
      lyt.plotOut();
    }

    // FREE MEMORY
    if ( OPTS::HIGHLIGHT ) {
      delete hLayer;
      hLayer = 0;
    }
    if ( OPTS::GC && OPTS::HASREF ) {
      delete gcLayer;
      gcLayer = 0;
    }

  }
  delete pRefLen;
  pRefLen = 0;
  return 0;
}

bool computeBAMCoverage(std::string& bamF, PrePlot& pPlot, unsigned int& refLength, std::vector<Cov>& layers) {
  // Opens a BAM FILE
  char *symlink( (char*) bamF.c_str());
  char *bamFilePath;
  bamFilePath = realpath(symlink, NULL);
  if (bamFilePath == NULL) {
    free(bamFilePath);
    std::cout << "ERROR: Cannot resolve filepath: " << bamF << std::endl;
    return false;
  }

  // Sets variables
  std::string refName( pPlot.getName() );
  std::string baiFilePath(std::string(bamFilePath)+".bai");
  seqan::CharString rName(refName);
  seqan::CharString bamFileName(bamFilePath);
  seqan::CharString baiFileName(baiFilePath);

  // Opens the BAI index.
  seqan::BamIndex<seqan::Bai> baiIndex;
  if ( !open(baiIndex, seqan::toCString(baiFileName)) ) {
      std::cout << "ERROR: Could not read BAI index file: " << baiFileName << "\n";
      return false;
  }

  // Open BamFileIn for reading.
  seqan::BamFileIn bamIn;
  if ( !open(bamIn, seqan::toCString(bamFileName)) ) {
      std::cout << "ERROR: Could not open BAM file: " << bamF << std::endl;
      return false;
  }

  // Reads header
  seqan::BamHeader header;
  seqan::readHeader(header, bamIn);

  // Translates from reference name to rID.
  int rID(0);
  if ( !seqan::getIdByName( rID, seqan::contigNamesCache( seqan::context(bamIn) ), rName ) ) {
    std::cout << "ERROR: Reference sequence named: " << refName << " not known" << std::endl;
    return false;
  }

  // SAM IS 1-INDEXED --> -1 at begin position
  // BED IS 0-INDEXED
  // Translate BEGIN and END arguments to number, 1-based (bed/commandline file) to 0-based (c++).
  // WARNING THIS IS NOT MOST EFFICIENT BUT IT IS THE SAFEST WAY TO DO IT
  // IDEA I WILL FIND A WAY TO BREAK WHEN PREPLOT END IS PASSED BY FAR
  unsigned int beginPos = 0;
  unsigned int endPos = refLength - 1; // REFERENCE LENGTH --> COVER WHOLE REFERENCE

  // Jump the BGZF stream to this position.
  bool hasAlignments(false);
  while ( !hasAlignments ) {
    if ( !seqan::jumpToRegion(bamIn, hasAlignments, rID, beginPos, endPos, baiIndex) ) {
        std::cout << "ERROR: Could not jump to " << beginPos << ":" << endPos << "\n";
        return false;
    }
    beginPos++; // GOES TO NEXT POSITION ON REFERENCE
  }
  // NOTE I WILL USE A SWITCH CASE FUNCTION AND RETURN INTEGERS 0: PASS 1: FAIL 2: SKIP ALIGNMENT -> DO NOT PLOT
  // NOTE DOING THIS WILL BE TIME CONSUMING SO NOT DONE YET
  if ( !hasAlignments ) {
    std::cout << "ERROR: no alignments found!" << std::endl;
    return false;
  }

  // WARNING NOW beginPos is at the first found alignment --> IT MAY **NOT** BE IN RANGE OF pPlot.getStart() : pPlot.getEnd()
  seqan::BamAlignmentRecord record;
  //seqan::BamFileOut out(bamIn, std::cout, seqan::Sam());

  // COVERAGE INFORMATION OF THE RECORD
  unsigned int bS( pPlot.getStart() );    // BED START
  unsigned int bE( pPlot.getEnd() );      // BED END
  // CHANGE BED END IN CASE IT WAS SET TO BIGGER THAN ACTUAL CONTIG SIZE
  if (bE > refLength) { bE = refLength - 1; } // SAM IS 0 INDEXED
  unsigned int lengthOfRegion(bE - bS);
  // CREATE A COVERAGE VECTOR AND FILLS IT WITH ZEROS
  std::vector<unsigned int> vCov(lengthOfRegion, 0);

  /*/ COMPUTING COVERAGE /*/
  do { // While not at end and current record ref == refName
    seqan::readRecord(record, bamIn);
    if (seqan::getContigName(record, bamIn) != rName) { // CHECK CONTIG NAME
      break;
    }
    if (record.qual < OPTS::PHRED) { // CHECK PHRED SCORE
      continue;
    }

    unsigned int aS(record.beginPos);
    unsigned int aL(seqan::length(record.seq));  // ALIGN LENGTH
    if ( aL == 1 ) {
      if (record.seq == "*") {
        continue;
      }
    }
    unsigned int aE(aS + aL);             // ALIGN END
    //std::cout << "record information: aS: " << aS << " aL: " << aL << " aE: " << aE << std::endl;

    // GET ALIGNMENT START POSITION AND SEQUENCE LENGTH
    if ( bS > aE || bE < aS ) {
      // MOST LIKELY CASE
      //                 bS          bE
      //                 |===========| BED
      // ================================================== REF
      //  |=========|                    |=================| ALIGNMENT
      //  aS        aE                   aS                aE
      continue;
    } else if ( bS >= aS && bS <= aE && bE >= aE) {
      // FIRST CASE
      //          bS              bS     bE            bE
      //          |======= or ====|======| or =========| BED
      // ================================================== REF
      //          |======================| ALIGNMENT
      //          aS                     aE
      // + 1 from vCov[0] to vCov[align End - bed start] --> Alignment end is comprised within bed (bE >= aE)
      for (int i(0); i < (aE-bS); i++) { vCov[i]++; }
    } else if ( bS < aS && bE >= aE ) {
      // SECOND CASE
      //          bS                     bE         bE
      //          |======================|==========| BED
      // ================================================== REF
      //               |=================| ALIGNMENT
      //               aS                aE
      // Alignment STARTS in bed and alignment ENDS before or at bed end --> vCov +1 between aS - bS and aE - bS
      for (int i(aS-bS); i < (aE-bS); i++) { vCov[i]++; }
    } else if ( bS < aS && bE < aE) {
      // THIRD CASE
      //   bedStart                          bedEnd
      //     |================================| BED
      // ================================================== REF
      //                            |=================| ALIGNMENT
      //                          AlnStart          AlnEnd
      // Alignment STARTS WITHIN BED and ENDS AFTER BED --> +1 from aS - bS to bE - bS
      for (int i(aS-bS); i < (bE-bS); i++) { vCov[i]++; }
    } else if ( bS > aS && bE < aE ) {
      // FOURTH CASE
      //                 bS          bE
      //                 |===========| BED
      // ================================================== REF
      //  |===============================================| ALIGNMENT
      //  aS                                              aE
      for (int i(0); i < vCov.size(); i++);
    }
  } while ( ! atEnd(bamIn) );

  // Sets the output file name
  std::string output(std::string( std::string(bamFilePath) + "_" + refName + "_" + toString(bS) + "-" + toString(bE) + ".cov") );
  writeCov(vCov, bS, output, refName);
  // Adds a coverage object to Cov
  double stats[4];
  computeStats(vCov, stats);
  Cov cC(vCov, bamF, refName, lengthOfRegion, stats[1], stats[0], stats[2], stats[3]);
  layers.push_back(cC);
  return true;
}

bool computeSAMCoverage(std::string& samF, PrePlot& pPlot, unsigned int& refLength, std::vector<Cov>& layers) {
  // Opens a SAM FILE
  char *symlink( (char*) samF.c_str());
  char *samFilePath;
  samFilePath = realpath(symlink, NULL);
  if (samFilePath == NULL) {
    free(samFilePath);
    std::cout << "ERROR: Cannot resolve filepath: " << samF << std::endl;
    return false;
  }

  // Sets variables
  std::string refName( pPlot.getName() );
  seqan::CharString rName(refName);
  seqan::CharString samFileName(samFilePath);

  // Open BamFileIn for reading.
  seqan::BamFileIn samIn;
  if ( !open(samIn, seqan::toCString(samFileName)) ) {
      std::cout << "ERROR: Could not open SAM file: " << samF << std::endl;
      return false;
  }

  // Reads header
  seqan::BamHeader header;
  seqan::readHeader(header, samIn);

  // Translates from reference name to rID.
  int rID(0);
  if ( !seqan::getIdByName( rID, seqan::contigNamesCache( seqan::context(samIn) ), rName ) ) {
    std::cout << "ERROR: Reference sequence named: " << refName << " not known" << std::endl;
    return false;
  }

  // SAM IS 1-INDEXED --> -1 at begin position
  // BED IS 0-INDEXED
  // Translate BEGIN and END arguments to number, 1-based (bed/commandline file) to 0-based (c++).
  // WARNING THIS IS NOT MOST EFFICIENT BUT IT IS THE SAFEST WAY TO DO IT
  // IDEA I WILL FIND A WAY TO BREAK WHEN PREPLOT END IS PASSED BY FAR
  unsigned int beginPos = 0;
  unsigned int endPos = refLength - 1; // REFERENCE LENGTH --> COVER WHOLE REFERENCE

  // WARNING NOW beginPos is at the first found alignment --> IT MAY **NOT** BE IN RANGE OF pPlot.getStart() : pPlot.getEnd()
  seqan::BamAlignmentRecord record;
  //seqan::BamFileOut out(bamIn, std::cout, seqan::Sam());

  // COVERAGE INFORMATION OF THE RECORD
  unsigned int bS( pPlot.getStart() );    // BED START
  unsigned int bE( pPlot.getEnd() );      // BED END
  // CHANGE BED END IN CASE IT WAS SET TO BIGGER THAN ACTUAL CONTIG SIZE
  if (bE > refLength) { bE = refLength - 1; } // SAM IS 0 INDEXED
  unsigned int lengthOfRegion(bE - bS);
  // CREATE A COVERAGE VECTOR AND FILLS IT WITH ZEROS
  std::vector<unsigned int> vCov(lengthOfRegion, 0);

  /*/ COMPUTING COVERAGE /*/
  do { // While not at end and current record ref == refName
    seqan::readRecord(record, samIn);
    if (seqan::getContigName(record, samIn) != rName) { // CHECK CURRENT REF NAME
      continue;
    }
    if (record.qual < OPTS::PHRED) { // CHECK PHRED SCORE
      continue;
    }
    unsigned int aS(record.beginPos);
    unsigned int aL(seqan::length(record.seq));  // ALIGN LENGTH
    unsigned int aE(aS + aL);             // ALIGN END
    //std::cout << "record information: aS: " << aS << " aL: " << aL << " aE: " << aE << std::endl;

    // GET ALIGNMENT START POSITION AND SEQUENCE LENGTH
    if ( bS > aE || bE < aS ) {
      // MOST LIKELY CASE
      //                 bS          bE
      //                 |===========| BED
      // ================================================== REF
      //  |=========|                    |=================| ALIGNMENT
      //  aS        aE                   aS                aE
      continue;
    } else if ( bS >= aS && bS <= aE && bE >= aE) {
      // FIRST CASE
      //          bS              bS     bE            bE
      //          |======= or ====|======| or =========| BED
      // ================================================== REF
      //          |======================| ALIGNMENT
      //          aS                     aE
      // + 1 from vCov[0] to vCov[align End - bed start] --> Alignment end is comprised within bed (bE >= aE)
      for (int i(0); i < (aE-bS); i++) { vCov[i]++; }
    } else if ( bS < aS && bE >= aE ) {
      // SECOND CASE
      //          bS                     bE         bE
      //          |======================|==========| BED
      // ================================================== REF
      //               |=================| ALIGNMENT
      //               aS                aE
      // Alignment STARTS in bed and alignment ENDS before or at bed end --> vCov +1 between aS - bS and aE - bS
      for (int i(aS-bS); i < (aE-bS); i++) { vCov[i]++; }
    } else if ( bS < aS && bE < aE) {
      // THIRD CASE
      //   bedStart                          bedEnd
      //     |================================| BED
      // ================================================== REF
      //                            |=================| ALIGNMENT
      //                          AlnStart          AlnEnd
      // Alignment STARTS WITHIN BED and ENDS AFTER BED --> +1 from aS - bS to bE - bS
      for (int i(aS-bS); i < (bE-bS); i++) { vCov[i]++; }
    } else if ( bS > aS && bE < aE ) {
      // FOURTH CASE
      //                 bS          bE
      //                 |===========| BED
      // ================================================== REF
      //  |===============================================| ALIGNMENT
      //  aS                                              aE
      for (int i(0); i < vCov.size(); i++);
    }
  } while ( ! atEnd(samIn) );

  // Sets the output file name to output coverage
  std::string output(std::string( std::string(samFilePath) + "_" + refName + "_" + toString(bS) + "-" + toString(bE) + ".cov") );
  writeCov(vCov, bS, output, refName);
  // Adds a coverage object to Cov
  double stats[4];
  computeStats(vCov, stats);
  Cov cC(vCov, samF, refName, lengthOfRegion, stats[1], stats[0], stats[2], stats[3]);
  layers.push_back(cC);
  return true;
}

void computeStats( const std::vector<unsigned int>& vCov, double (&stats)[4] ) {
  // Compute mean & stdv
  unsigned long sum(0);
  unsigned long length( vCov.size() );
  for (int i(0); i<length; i++) {
    sum+= vCov[i];
  } // RANGE BASED LOOP
  double mean(sum / length);
  long double accum(0.0);
  std::for_each( vCov.begin(), vCov.end(), [&](const double d) {
    accum += (d - mean) * (d - mean);
  } );
  double variance( accum / length );
  double stdv(sqrt(accum / length)); // Standard deviation of coverage
  auto result = std::minmax_element(vCov.begin(), vCov.end());

  stats[0] = vCov[result.first - vCov.begin()];   // MIN
  stats[1] = vCov[result.second - vCov.begin()];  // MAX
  stats[2] = mean;  // MEAN
  stats[3] = stdv;  // STDV

}

void printHelp() {
  std::cerr << "This is help!" << std::endl;
}

inline bool fileExist(const std::string& filepath) {
    return (access(filepath.c_str(), F_OK) != -1);
}

bool writeCov(std::vector<unsigned int>& vCov, unsigned int& bS, std::string& output, std::string& refName) {
  std::cout << "Writing coverage info to: " << output << std::endl;
  try {
    std::ofstream cov;
    cov.open(output, std::fstream::in | std::fstream::out);
    if ( !cov ) { // IF FILE DOE NOT EXISTS -> CREATE IT
      cov.open(output, std::fstream::in | std::fstream::out | std::fstream::trunc);
    }
    for (int i(0); i < vCov.size(); i++) {
      cov << refName << "\t" << bS+i+1 << "\t" << vCov[i] << std::endl; // +1 because sam files are 1-indexed and bed files are 0 indexed
    }
    cov.close();
    return true;
  } catch (std::exception const& e) {
    std::cerr << "WARNING: Could not output to .cov file!"
    << std::endl << "ERROR:" << e.what() << std::endl;
    return false;
  }
}

std::vector<std::string> split(std::string& input, const char * delimiters) {
    int inLen(input.length());    // Get current line length
    char cAr[inLen + 1];          // Creates a char array
    strcpy(cAr, input.c_str());   // Converts string to char array
    std::vector<std::string> out;
    char *tkn(strtok(cAr, delimiters));
    while (tkn != NULL) { out.push_back(std::string(tkn)); tkn = strtok(NULL, delimiters); }
    return out;
}

bool readHeader(std::string& bamF, std::map<std::string, unsigned int> * pRefLen, std::vector<std::string>& refNamesToPlot) {
  char *symlink( (char*) bamF.c_str());
  char *bamFilePath;
  bamFilePath = realpath(symlink, NULL);
  if (bamFilePath == NULL) {
    free(bamFilePath);
    std::cout << "ERROR: Cannot resolve filepath: " << bamF << std::endl;
    return false;
  }
  seqan::CharString bamFileName(bamFilePath);
  seqan::BamFileIn bamIn; // READS SAM OR BAM
  // Opens the file SAM or BAM format are readable
  if (! open(bamIn, seqan::toCString(bamFileName))) {
    std::cout << "ERROR: Could not open: " << bamFileName << std::endl;
    return false;
  }
  try {
    // Copy header.
    seqan::BamHeader header;
    seqan::readHeader(header, bamIn);
    typedef seqan::FormattedFileContext<seqan::BamFileIn, void>::Type TBamContext;
    TBamContext const & bamContext = seqan::context(bamIn);

    // Vector of ref in headers to check in next step if all are found :
    std::vector<std::string> allRefsInHeaders;
    // Iterates header references to set length in map pRefLen
    std::string currentRefName; // Inits a std::string for seqan conversion (seqan::move)
    std::map<std::string, unsigned int>::iterator it; // Iterator over a map
    for (unsigned i = 0; i < seqan::length(seqan::contigNames(bamContext)); i++) {
      seqan::move(currentRefName, seqan::contigNames(bamContext)[i]); // Converts CharString to std::string
      allRefsInHeaders.push_back(currentRefName);
      it = (*pRefLen).find(currentRefName); // Finds the reference name
      if ( it != (*pRefLen).end() ) { // If reference is found in header
        unsigned int length(seqan::contigLengths(bamContext)[i]); // Gets length from header (0 == NOT FOUND)
        // Since we have multiple bam files length must be the same for reference contigs!!!
        if ( (*pRefLen)[currentRefName] != 0 && (*pRefLen)[currentRefName] != length) {
          std::cout << "ERROR: Reference contigs length not consistent across input .bam files!" << std::endl;
          return false;
        } else if ( (*pRefLen)[currentRefName] == 0 ) {
          (*pRefLen)[currentRefName] = length; // Sets length to ref in map
        }
      }
      currentRefName.clear(); // Resets the string
    }
    // NOTE SHOULD FIND A WAY TO ITERATE ONLY ONCE BUT I HAVE NOT YET
    // Iterates refNamesToPlot to check they are found in header if not return false
    for (int i(0); i < refNamesToPlot.size(); i++) {
      if (std::find(allRefsInHeaders.begin(), allRefsInHeaders.end(), refNamesToPlot[i]) == allRefsInHeaders.end()) {
        std::cout << "ERROR: " << refNamesToPlot[i] << " was not found in: " << bamF << std::endl;
        return false;
      }
    }
  } catch (seqan::Exception const & e) {
    std::cout << "ERROR: " << e.what() << std::endl;
    return false;
  }
  return true;
}

bool extractSeq(std::string& refFile, std::string& refName, unsigned int& pS, unsigned int& pE, std::string& region) {
  char *symlink( (char*) refFile.c_str() );
  char *seqFilePath;
  seqFilePath = realpath(symlink, NULL);
  if (seqFilePath == NULL) {
    free(seqFilePath);
    std::cout << "ERROR: Cannot resolve filepath:" << refFile << std::endl;
    return false;
  }
  seqan::CharString refFilePath(refFile);
  seqan::SeqFileIn seqIn; // READS SEQ FILES (FA FQ GENBANK...)
  seqan::CharString id;   // Ctg id
  seqan::Dna5String seq;  // sequence data
  // Tries opening
  if (! open(seqIn, seqan::toCString(refFilePath))) {
    std::cout << "ERROR: Could not open: " << refFilePath << std::endl;
    return false;
  }
  // Tries to read a sequence!
  try {
    do {
      seqan::readRecord(id, seq, seqIn);
      if (id == refName) {
        unsigned int length(pE - pS);
        std::string totalSeq;
        seqan::assign(totalSeq, seq);
        region = totalSeq.substr(pS, length);
        return true;
      } else {
        continue;
      }
    } while ( !atEnd(seqIn) );
  } catch (seqan::Exception const & e) {
    std::cout << "ERROR: " << e.what() << std::endl;
    return false;
  }
  std::cout << "ERROR: Could not find reference contig in: " << refFile << std::endl;
  return false;
}

float computeGC(std::string& dna) {
  int GC(0);
  int AT(0);
  int N(0);
  for (int i(0); i<dna.length(); i++) {
    switch(dna[i]) {
      case 'A': AT++; break;
      case 'T': AT++; break;
      case 'U': AT++; break;
      case 'a': AT++; break;
      case 't': AT++; break;
      case 'G': GC++; break;
      case 'C': GC++; break;
      case 'g': GC++; break;
      case 'c': GC++; break;
      default: N++;
    }
  }
  return ( (float) GC / (GC+AT) );
}

bool readHighlight(std::vector<Highlight>& vH, std::string& bedF, unsigned int& pS, unsigned int& pE) {
  char *bedsymlink( (char*) bedF.c_str() );
  char *bedFilePath;
  bedFilePath = realpath(bedsymlink, NULL);
  if (bedFilePath == NULL) {
    free(bedFilePath);
    std::cout << "ERROR: Cannot resolve filepath:" << bedF << std::endl;
    return false;
  }
  try {
    // Open input bed file.
    seqan::BedFileIn bedIn(seqan::toCString(bedFilePath));
    seqan::BedRecord<seqan::Bed3> record;
    // Read record by record
    while (!atEnd(bedIn)) {
      seqan::readRecord(record, bedIn); // reads the next record
      std::string refName;  // Inits a name variable
      std::string geneName; // Inits a name variable
      std::string sData;    // String that will get back data from record
      seqan::move(refName, record.ref); // Converts CharString to std::string
      seqan::assign(sData, record.data); // Converts CharString to std::string
      Color cCol; // inits a color object

      // GETS COLOR OF HIGHLIGHT
      int stringLength(sData.length()); // Get current line length
  		char c_array[stringLength + 1];   // Creates a char array
  		strcpy(c_array, sData.c_str());   // Converts string to char array
      int pos(0);
      char *token(strtok(c_array, "\t"));
      while (token != NULL) {
        switch(pos) {
          case 0: geneName = std::string(token); break;
          case 1:
          {
            int colors[3];
            int i(0);
            char *tkn(strtok(token, ","));
            while (tkn != NULL) {
              colors[i] = strtol(tkn, NULL, 10);
              tkn = strtok(NULL, ",");
              i++;
            }
            cCol = Color(colors[0],colors[1],colors[2]);
          }
          break;
        }
        token = strtok(NULL, "\t");
        pos++;
      }
      //std::cout << "\tHIGHLIGHT: " << refName << "\t" << geneName << "\t" << cCol.r << "," << cCol.g << "," << cCol.b << std::endl;
      // CHECK IF HIGHLIGHT IS IN SPAN: MULTIPLE CASES
      // 1 CONVERT COORDINATES TO NEW REGION
      bool addHglt(true);
      unsigned int nHS(record.beginPos - pS);
      unsigned int nHE(record.endPos - pS);
      if (nHS < 0) {  // CASE START BEFORE SPAN
        nHS = 0;
      } else if (nHS > pE) { // CASE START AFTER SPAN
        addHglt = false;
      }

      if (nHE > pE - pS) { // IN CASE nHE > length of to plot region --> say nHE = end of region
        nHE = pE - pS;
      } else if (nHE < 0) { // IN CASE nHE is out of span
        addHglt = false;
      }

      if (addHglt) {
        // Create a highlight object
        Highlight hglt(refName, geneName, nHS, nHE, cCol);
        vH.push_back(hglt); // Pushes the highlight to vector reference
      }
    }
  } catch (seqan::Exception const & e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return false;
  }
  return true;
}
