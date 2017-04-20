#include <iostream>
#include <cstring>
#include <vector>
#include <numeric>
#include "../../tools/tinyxml2/tinyxml2.h"

using namespace std;

typedef uint16_t CountType;

struct SequenceElement {

	/** The id of the read. Equals the position in the argument list and in the AlignmentSettings::seqs vector (0-based). */
	CountType id;
	/** The mate number. 0 for barcodes, increasing for sequence reads in the given order (1-based). */
	CountType mate;
	/** The length of the respective read. */
	CountType length;

	/**
	 * Constructor of a SequenceElement NULL object.
	 * @author Tobias Loka
	 */
	SequenceElement () : id(0), mate(0), length(0) {};

	/**
	 * Constructor of a valid SequenceElement object.
	 * @param id The id of the read.
	 * @param m The mate number of the read (0 for barcodes, incrementing for sequence reads)
	 * @param l The length of the read
	 * @author Tobias Loka
	 */
	SequenceElement (CountType id, CountType m, CountType l): id(id), mate(m), length(l) {};

	/**
	 * Check whether the SequenceElement object is a barcode or not.
	 * @return true, if SequenceElement is a barcode. False if not.
	 * @author Tobias Loka
	 */
	bool isBarcode() { return (mate==0);}
};

struct AlignmentSettings {
  std::vector<uint16_t> lanes;
  std::vector<uint16_t> tiles;
  std::vector<SequenceElement> seqs;
  uint16_t mates = 0;
};



int main(int argc, char *argv[])
{
  AlignmentSettings settings;

  // Parse RunInfo.xml if present
  tinyxml2::XMLDocument doc;
  //if (!doc.LoadFile( settings.runInfo_fname ))
  if (!doc.LoadFile( "/home/schulzeja/private/RunInfo.xml" ))
    if (doc.FirstChildElement("RunInfo"))
      if (doc.FirstChildElement("RunInfo")->FirstChildElement("Run")) {
        if (doc.FirstChildElement("RunInfo")->FirstChildElement("Run")->FirstChildElement("Reads")) {
          tinyxml2::XMLElement * reads = doc.FirstChildElement("RunInfo")->FirstChildElement("Run")->FirstChildElement("Reads");
          if (reads->FirstChildElement("Read")) {
            tinyxml2::XMLElement * read = reads->FirstChildElement("Read");
            do {
              settings.seqs.push_back(SequenceElement(settings.seqs.size(), (*(read->Attribute("IsIndexedRead")) == 'N') ? ++settings.mates : 0, read->IntAttribute("NumCycles")));
            } while(read = read->NextSiblingElement());
          }
        }
        if (doc.FirstChildElement("RunInfo")->FirstChildElement("Run")->FirstChildElement("FlowcellLayout")) {
          tinyxml2::XMLElement * flowcellLayout = doc.FirstChildElement("RunInfo")->FirstChildElement("Run")->FirstChildElement("FlowcellLayout");

          std::vector<uint16_t> temp(flowcellLayout->IntAttribute("LaneCount"));
          std::iota(temp.begin(), temp.end(), 1);
          settings.lanes = temp;

          std::vector<uint16_t> temp2;
          for (uint16_t l = 1; l <= flowcellLayout->IntAttribute("SurfaceCount"); l++)
            for (uint16_t s = 1; s <= flowcellLayout->IntAttribute("SwathCount"); s++)
              for (uint16_t t = 1; t <= flowcellLayout->IntAttribute("TileCount"); t++)
                temp2.push_back( l*1000 + s*100 + t );
          settings.tiles = temp2;
        }
      }

  std::cout << "lanes" << std::endl;
  for (auto e:settings.lanes)
    std::cout << e << std::endl;
  std::cout << "tiles" << std::endl;
  for (auto e:settings.tiles)
    std::cout << e << std::endl;
  std::cout << "seqs" << std::endl;
  for (auto e:settings.seqs)
    std::cout << e.id << "\t" << e.mate << "\t" << e.length << std::endl;
  std::cout << "mates = " << settings.mates << std::endl;



    //tinyxml2::XMLDocument doc;
    //if (doc.LoadFile( "/home/schulzeja/private/RunInfo.xml" ))
        //return 1;

    //tinyxml2::XMLElement * reads = doc.FirstChildElement("RunInfo")->FirstChildElement("Run")->FirstChildElement("Reads");
    //{
    //std::cout << "ehem .." << std::endl;
    //tinyxml2::XMLElement * jak;
    //jak = doc.FirstChildElement("bla");
    //if (!jak)
      //break;

    //std::cout << "ehem 1.." << std::endl;
    //jak = jak->FirstChildElement("foen");
    //if (!jak)
      //break;

    //std::cout << "ehem 2.." << std::endl;
    //jak = jak->FirstChildElement("Reads");
    //if (!jak)
      //break;
    //}

    //int seqs[2] = {0, 0};
    //int bars[2] = {0, 0};
    //int number;

    //int seq1;
    //int seq2;
    //int bar1;
    //int bar2;
    //tinyxml2::XMLElement * read = reads->FirstChildElement("Read");

    //do {
        //number = read->IntAttribute("Number") - 1;
        //std::cout << "===============" << std::endl;
        //std::cout << read->Attribute("IsIndexedRead") << std::endl;
        //std::cout << std::strcmp(read->Attribute("IsIndexedRead"), "N") << std::endl;
        //std::cout << read->IntAttribute("NumCycles") << std::endl;
        //if (std::strcmp(read->Attribute("IsIndexedRead"), "N")==0)
            //seqs[number] = read->IntAttribute("NumCycles");
        //else
            //bars[number] = read->IntAttribute("NumCycles");
    //} while(read = read->NextSiblingElement());

    //std::cout << "seq1: " << seqs[1] << std::endl;
    //std::cout << "seq2: " << seqs[2] << std::endl;
    //std::cout << "bar1: " << bars[1] << std::endl;
    //std::cout << "bar2: " << bars[2] << std::endl;
    
    return 0;
}
