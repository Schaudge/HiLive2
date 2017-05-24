#include "../lib/headers.h"
#include "../lib/definitions.h"
#include "../lib/global_variables.h"
#include "../lib/kindex.h"
#include "../lib/alnstream.h"
#include "../lib/parallel.h"
#include "../lib/argument_parser.h"

int main() {

  seqan::CharString samFileName = "writtenBam.bam";
  seqan::BamFileOut samFileOut(seqan::toCString(samFileName));
  seqan::StringSet<seqan::CharString> referenceNames;
  //seqan::appendValue(referenceNames, "dummy");
  seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > referenceNamesCache(referenceNames);
  seqan::BamIOContext<seqan::StringSet<seqan::CharString> > bamIOContext(referenceNames, referenceNamesCache);
  seqan::contigNames(bamIOContext) = std::vector<std::string>{"gi|568815576|ref|NC_000022.11|", "gi|568815576|ref|NC_000022.11|"};
  seqan::contigLengths(bamIOContext) = std::vector<unsigned>{568815576, 568815576};
  samFileOut.context = bamIOContext;


  /////////////////
  // set SAM header
  seqan::BamHeaderRecord headerRecord;
  
  // @HD header.
  seqan::clear(headerRecord);
  headerRecord.type = seqan::BAM_HEADER_FIRST;
  seqan::resize(headerRecord.tags, 2);
  headerRecord.tags[0].i1 = "VN";
  headerRecord.tags[0].i2 = "1.5";
  headerRecord.tags[1].i1 = "GO";
  headerRecord.tags[1].i2 = "query";
  seqan::writeHeader(samFileOut, headerRecord);

  // @SQ header.
  //std::stringstream ss;
  //ss.str(std::string()); // clear string stream
  //seqan::clear(headerRecord);
  //headerRecord.type = seqan::BAM_HEADER_REFERENCE;
  //seqan::resize(headerRecord.tags, 2);
  //headerRecord.tags[0].i1 = "SN";
  //headerRecord.tags[0].i2 = "gi|568815576|ref|NC_000022.11|";

  //headerRecord.tags[1].i1 = "LN";
  //ss << 50818468;
  //headerRecord.tags[1].i2 = ss.str();

  //seqan::writeHeader(samFileOut, headerRecord);
  std::cout << "Checkpoint A" << std::endl;
  return 0;


  //std::string filename = "/scratch/schulzeja/workDir/benchlive/runs/Y2017M05D10_14.22.29_jak_runXmlParsing_pairedEnd/samFiles/L001/s_1_1101.1.bam";
  //std::string filename = "/scratch/schulzeja/workDir/benchlive/runs/Y2017M05D10_14.03.05_jak_runXmlParsing_pairedEnd/samFiles/L001/s_1_1101.1.sam";
  
  // without PG header
  std::string filename = "/scratch/schulzeja/workDir/benchlive/runs/Y2017M05D10_15.42.33_jak_runXmlParsing_pairedEnd/samFiles/L001/s_1_1101.1.bam";

  seqan::BamAlignmentRecord record1;
  seqan::BamAlignmentRecord record2;
  seqan::BamHeader header1;
  seqan::BamHeader header2;
  //seqan::BamHeader header3;

  std::cout << "Reading " << filename << " ..." << std::endl; // Checkpoint for debugging
  seqan::BamFileIn bamFileIn(filename.c_str());





  //seqan::BamHeader header;
  //seqan::readHeader(header, bamFileIn);

  //typedef seqan::FormattedFileContext<seqan::BamFileIn, void>::Type TBamContext;

  //TBamContext const & bamContext = seqan::context(bamFileIn);

  //std::cout << seqan::length(seqan::contigNames(bamContext)) << std::endl;
  //for (unsigned i = 0; i < seqan::length(seqan::contigNames(bamContext)); ++i)
      //std::cout << seqan::contigNames(bamContext)[i] << '\t' << seqan::contigLengths(bamContext)[i] << '\n';





  std::cout << "Checkpoint 8.1" << std::endl;
  seqan::readHeader(header1, bamFileIn);
  seqan::readHeader(header2, bamFileIn);
  //seqan::readHeader(header3, bamFileIn);

  std::cout << "Checkpoint 8.2" << std::endl;
  seqan::readRecord(record1, bamFileIn);
  seqan::readRecord(record2, bamFileIn);





  //seqan::writeRecord(samFileOut, record1);





  std::cout << "Checkpoint 8.3" << std::endl;
  seqan::BamFileOut bamFileOut1(seqan::context(bamFileIn), "dummySam1.sam");
  seqan::writeHeader(bamFileOut1, header1);
  seqan::BamFileOut bamFileOut2(seqan::context(bamFileIn), "dummySam2.sam");
  seqan::writeHeader(bamFileOut2, header2);
  //seqan::BamFileOut bamFileOut3(seqan::context(bamFileIn), "dummySam3.sam");
  //seqan::writeHeader(bamFileOut3, header3);

  seqan::BamFileOut bamFileOut4(seqan::context(bamFileIn), "dummySam4.sam");
  seqan::writeRecord(bamFileOut4, record1);
  seqan::BamFileOut bamFileOut5(seqan::context(bamFileIn), "dummySam5.sam");
  seqan::writeRecord(bamFileOut5, record2);
}
