#include "bed.h"

// Classes functions
std::vector<PrePlot> readBed(const std::string& bedFilePath, std::map<std::string, unsigned int> * pRefLen) {
  // Create a preplot object vector
  std::vector<PrePlot> out;
  try {
    // Open input bed file.
    seqan::BedFileIn bedIn(seqan::toCString(bedFilePath));
    // Create output stream to bed file.
    seqan::BedFileOut bedOut(std::cout, seqan::Bed());

    // Copy the file record by record.
    seqan::BedRecord<seqan::Bed3> record;

    // Read record by record
    while (!atEnd(bedIn)) {
      seqan::readRecord(record, bedIn); // reads the next record
      std::string name; // Inits a name variable
      seqan::move(name, record.ref); // Converts CharString to std::string
      PrePlot pPlt(name, (unsigned int) record.beginPos, (unsigned int) record.endPos); // Create a preplot object
      (*pRefLen)[name] = (unsigned int) 0; // Inits a key name (reference) to 0 length (will be found in .bai or in header context)
      out.push_back(pPlt); // Fills the vector to return
    }
  }
  catch (seqan::Exception const & e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
  }

  return out;
}
