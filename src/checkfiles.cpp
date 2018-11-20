#include "checkfiles.h"

// FUNCTIONS
bool checkBAM(std::string& bamF) {
  // Try OPENING BAM FILE
  char *symlink( (char*) bamF.c_str());
  char *bamFilePath;
  bamFilePath = realpath(symlink, NULL);
  if (bamFilePath == NULL) {
    free(bamFilePath);
    std::cout << "ERROR: Cannot resolve filepath: " << bamF << std::endl;
    return false;
  }

  // CHECK IF FILES ARE INDEXED (BAI FILES ARE FOUND)
  std::string baiFilePath(std::string(bamFilePath)+".bai");
  if ( !fileExist(baiFilePath) ) {
    std::cout << "ERROR: Cannot find index file: " << baiFilePath << "!" << std::endl;
    return false;
  }
  // CHECK IF FILES ARE SORTED
  seqan::CharString bamFileName(bamFilePath);
  seqan::BamFileIn bamIn; // READS SAM OR BAM
  // Opens the file SAM or BAM format are readable
  if (! open(bamIn, seqan::toCString(bamFileName))) {
    std::cout << "ERROR: Could not open: " << bamFileName << std::endl;
    return false;
  }
  try {
    // Copy header.
    seqan::BamHeader header; // Class- BamHeader is a String<BamHeaderRecord>
    seqan::readHeader(header, bamIn);
    seqan::BamHeaderRecord record;
    bool isSorted(false);
    bool isSortedByCoordinates(false);
    std::string msg;
    seqan::CharString sortedTag("SO");
    seqan::CharString value;
    seqan::CharString cCtgName;
    for (int i(0); i < seqan::length(header); i++) {
      record = header[i]; // set that record =
      if ( seqan::getTagValue(value, sortedTag, record) ) {
        isSorted = true;
        if ( value == seqan::CharString("coordinate") ) {
          isSortedByCoordinates = true;
          break; // Break the loop as not needed any more!
        }
      }
    }
    if ( !isSorted ) {
      std::cout << "ERROR: BAM file is not sorted!" << std::endl;
      return false;
    } else if ( !isSortedByCoordinates ) {
      std::cout << "ERROR: BAM file is not sorted by coordinates!" << std::endl;
      return false;
    }
  } catch (seqan::Exception const & e) {
    std::cout << "ERROR: " << e.what() << std::endl;
    return false;
  }
  std::cout << "\tBAM: " << bamF << std::endl;
  return true; // All bams exist, are indexed and are sorted by coordinates
}

bool checkSAM(std::string& samF) {
  // Try OPENING SAM FILE
  char *symlink( (char*) samF.c_str());
  char *samFilePath;
  samFilePath = realpath(symlink, NULL);
  if (samFilePath == NULL) {
    free(samFilePath);
    std::cout << "ERROR: Cannot resolve filepath: " << samF << std::endl;
    return false;
  }

  // CHECK IF FILES ARE SORTED
  seqan::CharString samFileName(samFilePath);
  seqan::BamFileIn samIn; // READS SAM OR BAM
  // Opens the file SAM or BAM format are readable
  if (! open(samIn, seqan::toCString(samFileName))) {
    std::cout << "ERROR: Could not open: " << samFileName << std::endl;
    return false;
  }
  try {
    // Copy header.
    seqan::BamHeader header; // Class- BamHeader is a String<BamHeaderRecord>
    seqan::readHeader(header, samIn);
    seqan::BamHeaderRecord record;
    bool isSorted(false);
    bool isSortedByCoordinates(false);
    std::string msg;
    seqan::CharString sortedTag("SO");
    seqan::CharString value;
    seqan::CharString cCtgName;
    for (int i(0); i < seqan::length(header); i++) {
      record = header[i]; // set that record =
      if ( seqan::getTagValue(value, sortedTag, record) ) {
        isSorted = true;
        if ( value == seqan::CharString("coordinate") ) {
          isSortedByCoordinates = true;
          break; // Break the loop as not needed any more!
        }
      }
    }
    if ( !isSorted ) {
      std::cout << "ERROR: SAM file is not sorted!" << std::endl;
      return false;
    } else if ( !isSortedByCoordinates ) {
      std::cout << "ERROR: SAM file is not sorted by coordinates!" << std::endl;
      return false;
    }
  } catch (seqan::Exception const & e) {
    std::cout << "ERROR: " << e.what() << std::endl;
    return false;
  }
  std::cout << "\tSAM: " << samF << std::endl;
  return true; // All bams exist, are indexed and are sorted by coordinates
}

bool checkREF(std::string& ref) {
  // CHECK IF FILE EXISTS
  // Try OPENING SEQ FILE
  // Creates a seqan string for storing .bam filename
  char *symlink( (char*) ref.c_str() );
  char *seqFilePath;
  seqFilePath = realpath(symlink, NULL);
  if (seqFilePath == NULL) {
    free(seqFilePath);
    std::cout << "ERROR: Cannot resolve filepath:" << ref << std::endl;
    return false;
  }
  seqan::CharString refFilePath(ref);
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
    seqan::readRecord(id, seq, seqIn);
  } catch (seqan::Exception const & e) {
    std::cout << "ERROR: " << e.what() << std::endl;
    return false;
  }
  std::cout << "\tREF: " << ref << std::endl;
  return true; // ALL PASSED
}

bool checkBED(std::string& bed) {
  char *bedsymlink( (char*) bed.c_str() );
  char *bedFilePath;
  bedFilePath = realpath(bedsymlink, NULL);
  if (bedFilePath == NULL) {
    free(bedFilePath);
    std::cout << "ERROR: Cannot resolve filepath:" << bed << std::endl;
    return false;
  }
  try {
    // Open input bed file.
    seqan::BedFileIn bedIn(seqan::toCString(bedFilePath));
    seqan::BedRecord<seqan::Bed3> record;
    // Read record by record
    seqan::readRecord(record, bedIn); // reads the next record
  }
  catch (seqan::Exception const & e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return false;
  }
  std::cout << "\tBED: " << bed << std::endl;
  return true;
}

inline bool fileExist(const std::string& filepath) {
    return (access(filepath.c_str(), F_OK) != -1);
}
