#include "argument_parser.h"
namespace po = boost::program_options;


int parseCommandLineArguments(std::string license, int argc, char const ** argv)
{
    po::options_description general("General");
    general.add_options()
        ("help,h", "Print this help message and exit")
        ("license", "Print licensing information and exit");

    po::options_description parameters("Parameters");
    parameters.add_options()
        ("BC_DIR", po::value<std::string>()->required(), "Illumina BaseCalls directory")
        ("INDEX", po::value<std::string>()->required(), "Path to k-mer index")
        ("CYCLES", po::value<CountType>()->required(), "Number of cycles")
        ("OUTDIR", po::value<std::string>(), "Directory to store sam files in [Default: temporary or BaseCalls directory");

    po::options_description io_settings("IO settings");
    io_settings.add_options()
        ("temp", po::value<std::string>()->default_value(""), "Temporary directory for the alignment files [Default: use BaseCalls directory]")
        ("bam,B", po::bool_switch()->default_value(false), "Create BAM files instead of SAM files [Default: false]")
        ("extended-cigar", po::bool_switch()->default_value(false), "Activate extended CIGAR format (= and X instead of only M) in output files [Default: false]")
        ("keep-files,k", po::bool_switch()->default_value(false), "Keep intermediate alignment files [Default: false]")
        ("lanes,l", po::value< std::vector<uint16_t> >()->multitoken()->composing(), "Select lane [Default: all lanes]")
        ("tiles,t", po::value< std::vector<uint16_t> >()->multitoken()->composing(), "Select tile numbers [Default: all tiles]")
        ("reads,r", po::value< std::vector<std::string> >()->multitoken()->composing(), "Enumerate read lengths and type. Example: -r 101R 8B 8B 101R equals paired-end sequencing with 2x101bp reads and 2x8bp barcodes. Overwrites information of runInfo.xml. [Default: single end reads without barcodes]")
        ("runInfoPath", po::value<std::string>(), "Path to runInfo.xml for parsing read and index lengths [Default: BC_DIR/../../RunInfo.xml]")
        ("barcodes,b", po::value< std::vector<std::string> >()->multitoken()->composing(), "Enumerate barcodes (must have same length) for demultiplexing, e.g. -b AGGATC -b CCCTTT [Default: no demultiplexing]")
    	("barcode-errors,E", po::value< std::vector<uint16_t> >()->multitoken()->composing(), "Enumerate the number of tolerated errors (only SNPs) for each barcode fragment, e.g. -E 2 2 [Default: 1 per fragment]")
		("keep-all-barcodes", po::bool_switch()->default_value(false), "Align and output all barcodes [Default: false]");

    po::options_description alignment("Alignment settings");
    alignment.add_options()
        ("min-errors,e", po::value<CountType>()->default_value(2), "Number of errors tolerated in read alignment [Default: 2]")
        ("all-best-hit,H", po::bool_switch()->default_value(false), "Report all of the best alignments for each read")
        ("any-best-hit", po::bool_switch(), "Report one of the best alignments for each read (default)")
        ("all-best-n-scores,N", po::value<CountType>(), "Report all alignments of the N best alignment scores for each read")
        ("all-hits,A", po::bool_switch()->default_value(false), "Report all valid alignments for each read")
        ("disable-ohw-filter", po::bool_switch()->default_value(true), "disable the One-Hit Wonder filter [Default: false]")
        ("start-ohw", po::value<CountType>()->default_value(globalAlignmentSettings.get_kmer_weight()+5), "First cycle to apply One-Hit Wonder filter [Default: K+5]")
        ("window,w", po::value<DiffType>()->default_value(5), "Set the window size to search for alignment extension, i.e. maximum total insertion/deletion size [Default: 5]")
        ("min-quality", po::value<CountType>()->default_value(1), "Minimum allowed basecall quality [Default: 1]");

    po::options_description technical("Technical settings");
    technical.add_options()
        ("block-size", po::value<uint64_t>()->default_value(64*1024*1024), "Block size for the alignment input/output stream in Bytes. Use -K or -M to specify in Kilobytes or Megabytes")
        (",K", po::bool_switch()->default_value(false), "Interpret the block-size argument as Kilobytes instead of Bytes")
        (",M", po::bool_switch()->default_value(false), "Interpret the block-size argument as Megabytes instead of Bytes")
        ("compression,c", po::value<uint8_t>()->default_value(2), "Compress alignment files. 0: no compression (default) 1: Deflate (smaller) 2: LZ4 (faster)")
        ("num-threads,n", po::value<int>(), "Number of threads to spawn [Default: all available]");

    po::options_description cmdline_options;
    cmdline_options.add(general).add(parameters).add(io_settings).add(alignment).add(technical);

    po::options_description visible_options;
    visible_options.add(general).add(io_settings).add(alignment).add(technical);

    std::stringstream help_message;
    help_message << "HiLive v"<< HiLive_VERSION_MAJOR << "." << HiLive_VERSION_MINOR << " - Realtime Alignment of Illumina Reads" << std::endl;
    help_message << "Copyright (c) 2015, Martin S. Lindner" << std::endl;
    help_message << "HiLive is open-source software. Check with --license for details." << std::endl << std::endl;
    help_message << "Usage: " << std::string(argv[0]) << " [options] BC_DIR INDEX CYCLES OUTDIR" << std::endl;
    help_message << "  BC_DIR       Illumina BaseCalls directory of the sequencing run to analyze" << std::endl;
    help_message << "  INDEX        Path to k-mer index file (*.kix)" << std::endl;
    help_message << "  CYCLES       Total number of cycles for read 1" << std::endl;
    help_message << "  OUTDIR       Directory to store sam files in" << std::endl;
    help_message << visible_options << std::endl;

    po::positional_options_description p;
    p.add("BC_DIR", 1);
    p.add("INDEX", 1);
    p.add("CYCLES", 1);
    p.add("OUTDIR", 1);

    
    // parse the arguments

    po::variables_map vm;
    try {
        // parse arguments
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
        // first check if -h or --help was called
        if (vm.count("help")) {
            std::cout << help_message.str();
            return 1;
        }
        // first check if --license was called
        if (vm.count("license")) {
            std::cout << license << std::endl;
            return 1;
        }

        // then check arguments
        po::notify(vm);  
    }
    catch ( boost::program_options::required_option& e ) {
        std::cerr << "Missing Parameter: " << e.what() << std::endl << std::endl;
        std::cout << help_message.str();
        return -1;  
    }
    catch( boost::program_options::error& e) { 
        std::cerr << "Error while parsing command line options: " << e.what() << std::endl << std::endl; 
        std::cout << help_message.str();
        return -1;  
    } 


    globalAlignmentSettings.set_root(vm["BC_DIR"].as<std::string>());
    globalAlignmentSettings.set_index_fname(vm["INDEX"].as<std::string>());
    globalAlignmentSettings.set_cycles(vm["CYCLES"].as<CountType>());
    globalAlignmentSettings.set_temp_dir(vm["temp"].as<std::string>());
    globalAlignmentSettings.set_write_bam(vm["bam"].as<bool>());
    globalAlignmentSettings.set_extended_cigar(vm["extended-cigar"].as<bool>());
    globalAlignmentSettings.set_keep_aln_files(vm["keep-files"].as<bool>());
    globalAlignmentSettings.set_keep_all_barcodes(vm["keep-all-barcodes"].as<bool>());
    globalAlignmentSettings.set_min_errors(vm["min-errors"].as<CountType>());

    globalAlignmentSettings.set_discard_ohw(vm["disable-ohw-filter"].as<bool>());
    globalAlignmentSettings.set_start_ohw(vm["start-ohw"].as<CountType>());
    globalAlignmentSettings.set_window(vm["window"].as<DiffType>());
    globalAlignmentSettings.set_min_qual(vm["min-quality"].as<CountType>());
    globalAlignmentSettings.set_compression_format(vm["compression"].as<uint8_t>());


    if (vm.count("OUTDIR"))
        globalAlignmentSettings.set_out_dir(boost::filesystem::path(vm["OUTDIR"].as<std::string>()));
    else {
        if (globalAlignmentSettings.get_temp_dir() == "") 
            globalAlignmentSettings.set_out_dir(boost::filesystem::path(globalAlignmentSettings.get_root()));
        else
            globalAlignmentSettings.set_out_dir(boost::filesystem::path(globalAlignmentSettings.get_temp_dir()));
    }


    // Parse RunInfo.xml if present
    if (vm.count("runInfoPath")) {
      globalAlignmentSettings.set_runInfo_fname(vm["runInfoPath"].as<std::string>());
      std::ifstream inputStream(globalAlignmentSettings.get_runInfo_fname());

      using boost::property_tree::ptree;
      ptree tree;

      read_xml(inputStream, tree);
      if (!tree.empty() && tree.count("RunInfo")!=0) {
        ptree ptree_RunInfo = tree.get_child("RunInfo");

        if (ptree_RunInfo.count("Run")!=0) {
          ptree ptree_Run = ptree_RunInfo.get_child("Run");

          if (ptree_Run.count("Reads")!=0) {
            ptree ptree_Reads = ptree_Run.get_child("Reads");
            unsigned lenSum = 0;
            std::vector<SequenceElement> tempSeqs;
            unsigned mates = 0; 
            for (const auto &read : ptree_Reads) {
              tempSeqs.push_back(SequenceElement(tempSeqs.size(), (read.second.get<std::string>("<xmlattr>.IsIndexedRead") == "N") ? ++mates : 0, read.second.get<unsigned>("<xmlattr>.NumCycles")));
              lenSum += read.second.get<unsigned>("<xmlattr>.NumCycles");
            }
            globalAlignmentSettings.set_seqs(tempSeqs);
            globalAlignmentSettings.set_mates(mates);
            if ( lenSum!=globalAlignmentSettings.get_cycles() ) {
              std::cerr << "Sum of parsed read lengths does not equal the given number of cycles." << std::endl;
              return false;
            }

            if (ptree_Run.count("FlowcellLayout")!=0) {
            ptree ptree_FlowcellLayout = ptree_Run.get_child("FlowcellLayout");

              std::vector<uint16_t> temp(ptree_FlowcellLayout.get<unsigned>("<xmlattr>.LaneCount"));
              std::iota(temp.begin(), temp.end(), 1);
              globalAlignmentSettings.set_lanes(temp);
              std::vector<uint16_t> temp2;
              for (uint16_t l = 1; l <= ptree_FlowcellLayout.get<unsigned>("<xmlattr>.SurfaceCount"); l++)
                for (uint16_t s = 1; s <= ptree_FlowcellLayout.get<unsigned>("<xmlattr>.SwathCount"); s++)
                  for (uint16_t t = 1; t <= ptree_FlowcellLayout.get<unsigned>("<xmlattr>.TileCount"); t++)
                    temp2.push_back( l*1000 + s*100 + t );
              globalAlignmentSettings.set_tiles(temp2);
            }
          }
        }
      }
    }
    

    // Parse read lengths and types.
    if ( vm.count("reads") && !parseReadsArgument(vm["reads"].as< std::vector<std::string> >()) )
    	return -1;
    else if ( !vm.count("reads") && globalAlignmentSettings.get_seqs().size() == 0 ) {
    	globalAlignmentSettings.set_seqs(std::vector<SequenceElement> {SequenceElement(0,1,globalAlignmentSettings.get_cycles())});
    	globalAlignmentSettings.set_mates(1);
    }

    // Parse lanes
    if (vm.count("lanes"))
        globalAlignmentSettings.set_lanes(vm["lanes"].as< std::vector<uint16_t> >());
    else
      if (globalAlignmentSettings.get_lanes().size() == 0) {
        std::vector<uint16_t> tempLanes = all_lanes();
        std::sort( tempLanes.begin(), tempLanes.end() );
        tempLanes.erase( std::unique( tempLanes.begin(), tempLanes.end() ), tempLanes.end() );
        globalAlignmentSettings.set_lanes(tempLanes);
      }

    // Parse tiles
    if (vm.count("tiles"))
        globalAlignmentSettings.set_tiles(vm["tiles"].as< std::vector<uint16_t> >());
    else
      if (globalAlignmentSettings.get_tiles().size() == 0) {
        std::vector<uint16_t> tempTiles = all_tiles();
        std::sort( tempTiles.begin(), tempTiles.end() );
        tempTiles.erase( std::unique( tempTiles.begin(), tempTiles.end() ), tempTiles.end() );
        globalAlignmentSettings.set_tiles(tempTiles);
      }


    if (vm.count("barcodes")) {
        if( !parseBarcodeArgument(vm["barcodes"].as< std::vector<std::string> >()) ) {
        	std::cerr << "Parsing error: Invalid barcode(s) detected. Please ensure that you used \'-\' "
        			"as duplex delimiter and that all barcodes have the correct length. Only use A,C,G and T as bases!" << std::endl;
        	return -1;
        }
    }

    if ( globalAlignmentSettings.get_barcodeVector().size() != 0 ) {
    	if ( vm.count("barcode-errors") ) {
    		globalAlignmentSettings.set_barcode_errors(vm["barcode-errors"].as< std::vector<uint16_t> >());
    		if ( globalAlignmentSettings.get_barcodeVector()[0].size() != globalAlignmentSettings.get_barcode_errors().size() ) {
    			std::cerr << "Parsing error: Number of barcode errors does not equal the number of barcodes." << std::endl;
    			return -1;
    		}
    	} else {
          std::vector<uint16_t> temp;
        	for ( uint16_t i = 0; i < globalAlignmentSettings.get_barcodeVector()[0].size(); i++ )
        		temp.push_back(1);
          globalAlignmentSettings.set_barcode_errors(temp);
    	}
    }

    if (vm["all-hits"].as<bool>()) {
        // all hits: disable other modes
        globalAlignmentSettings.set_any_best_hit_mode(false);
        globalAlignmentSettings.set_all_best_hit_mode(false);
        globalAlignmentSettings.set_all_best_n_scores_mode(false);
    } else if (vm["all-best-hit"].as<bool>()) {
        // enable all-best-hit mode and disable others
        globalAlignmentSettings.set_any_best_hit_mode(false);
        globalAlignmentSettings.set_all_best_hit_mode(true);
        globalAlignmentSettings.set_all_best_n_scores_mode(false);
    } else if (vm.count("all-best-n-scores")) {
        // enable any-best-n mode and get parameter
        globalAlignmentSettings.set_any_best_hit_mode(false);
        globalAlignmentSettings.set_all_best_hit_mode(false);
        globalAlignmentSettings.set_all_best_n_scores_mode(true);
        globalAlignmentSettings.set_best_n(vm["all-best-n-scores"].as<CountType>());
    } else { // the default behaviour
        // enable any-best-hit mode and disable others
        globalAlignmentSettings.set_any_best_hit_mode(true);
        globalAlignmentSettings.set_all_best_hit_mode(false);
        globalAlignmentSettings.set_all_best_n_scores_mode(false);
    }
        
    if (vm["-M"].as<bool>())
        globalAlignmentSettings.set_block_size(vm["block-size"].as<uint64_t>()*1024*1024);
    else if (vm["-K"].as<bool>())
        globalAlignmentSettings.set_block_size(vm["block-size"].as<uint64_t>()*1024);
    else
        globalAlignmentSettings.set_block_size(vm["block-size"].as<uint64_t>());

    if (vm.count("num-threads")) 
        globalAlignmentSettings.set_num_threads(vm["num-threads"].as<int>());
    else { // the default case, meaning as much as physically useful
        uint32_t n_cpu = std::thread::hardware_concurrency();
        if (n_cpu > 1)
            globalAlignmentSettings.set_num_threads(std::min( n_cpu, uint32_t( globalAlignmentSettings.get_lanes().size() * globalAlignmentSettings.get_tiles().size() ) ) );
        else
            globalAlignmentSettings.set_num_threads(1);
    }


    // check paths and file names
    if (!file_exists(globalAlignmentSettings.get_index_fname())){
        std::cerr << "Input error: Could not find k-mer index file " << globalAlignmentSettings.get_index_fname() << std::endl;
        return -1;
    }

    std::size_t found = globalAlignmentSettings.get_root().find("BaseCalls");
    if (!(found != std::string::npos && found >= globalAlignmentSettings.get_root().size()-10)) {
        std::cerr << "Warning: BaseCalls directory seems to be invalid: " << globalAlignmentSettings.get_root() << std::endl;
    } 

    if (!is_directory(globalAlignmentSettings.get_root())){
        std::cerr << "Input error: Could not find BaseCalls directory " << globalAlignmentSettings.get_root() << std::endl;
        return -1;
    }

    for ( uint16_t ln : globalAlignmentSettings.get_lanes() ) {
        std::string ln_dir = globalAlignmentSettings.get_root();
        if ( ln < 10 )
            ln_dir += "/L00";
        else if ( ln < 100 )
            ln_dir += "/L0";
        else
            ln_dir += "/L";
        ln_dir += std::to_string(ln);
        if (!is_directory(ln_dir)){
            std::cerr << "Input error: Could not find location of Lane " << ln << ": " << ln_dir << std::endl;
            return -1;
        }
    }


    // Compute maximal consecutive gaps in gap pattern
    CountType current_consecutive_gaps = 0;
    CountType last_gap = 0;
    CountType temp_max_consecutive_gaps = 0;

    for ( unsigned el : globalAlignmentSettings.get_kmer_gaps() ) {

    	// init first gap
    	if ( last_gap == 0 ) {
    		current_consecutive_gaps = 1;
    		last_gap = el;
    		continue;
    	}

    	// handle consecutive gaps
    	else if ( el == unsigned( last_gap + 1 ) ){
    		current_consecutive_gaps += 1;
    		last_gap = el;
    	}

    	// handle end of gap region
    	else {
    		temp_max_consecutive_gaps = std::max ( temp_max_consecutive_gaps, current_consecutive_gaps );
    		current_consecutive_gaps = 1;
    		last_gap = el;
    	}

    }
    globalAlignmentSettings.set_max_consecutive_gaps(std::max ( temp_max_consecutive_gaps, current_consecutive_gaps ) );


    // Report the basic settings
    std::cout << "Running HiLive with       " << globalAlignmentSettings.get_num_threads() << " thread(s)." << std::endl;
    std::cout << "BaseCalls directory:      " << globalAlignmentSettings.get_root() << std::endl;
    if (globalAlignmentSettings.get_temp_dir() != "") {
        std::cout << "Temporary directory:      " << globalAlignmentSettings.get_temp_dir() << std::endl;
    }
    if (!globalAlignmentSettings.get_write_bam())
        std::cout << "SAM output directory:     " << globalAlignmentSettings.get_out_dir() << std::endl;
    else
        std::cout << "BAM output directory:     " << globalAlignmentSettings.get_out_dir() << std::endl;
    std::cout << "Lanes:                    ";
    for ( uint16_t ln : globalAlignmentSettings.get_lanes() )
        std::cout << ln << " ";
    std::cout << std::endl;
    std::cout << "K-mer index:              " << globalAlignmentSettings.get_index_fname() << std::endl;
    std::cout << "Read lengths:             ";
    std::string barcode_suffix;
    for ( uint16_t read = 0; read != globalAlignmentSettings.get_seqs().size(); read ++) {
    	std::cout << globalAlignmentSettings.getSeqById(read).length;
    	barcode_suffix = globalAlignmentSettings.getSeqById(read).isBarcode() ? "B" : "R";
    	std::cout << barcode_suffix << " ";
    }
    std::cout << std::endl;
    std::cout << "Mapping error:            " << globalAlignmentSettings.get_min_errors() << std::endl;
    if (globalAlignmentSettings.get_any_best_hit_mode()) 
        std::cout << "Mapping mode:             Any-Best-Hit-Mode" << std::endl;
    else if (globalAlignmentSettings.get_all_best_hit_mode()) 
        std::cout << "Mapping mode:             All-Best-Hit-Mode" << std::endl;
    else if (globalAlignmentSettings.get_all_best_n_scores_mode()) 
        std::cout << "Mapping mode:             All-Best-N-Scores-Mode with N=" << globalAlignmentSettings.get_best_n() << std::endl;
    else
        std::cout << "Mapping mode:             All-Hits-Mode" << std::endl;
    std::cout << std::endl;

    return 0;
}


bool parseReadsArgument(std::vector< std::string > readsArg){

	CountType lenSum = 0;
	CountType length = 0;
	std::string length_string = "";
	char type;
	unsigned mates = 0;
  std::vector<SequenceElement> temp;

	for ( auto read = readsArg.begin(); read != readsArg.end(); ++read ) {

		length_string = (*read).substr(0,(*read).length()-1);
		type = (*(*read).rbegin());

		try{
		std::stringstream( length_string ) >> length;
		} catch( std::bad_cast & ex ){
			std::cerr << "Error while casting length " << length_string << " to type uint16_t" << std::endl;
		}

		if ( type!='B' && type!='R' ) {
			std::cerr << "\'" << type << "\'" << " is no valid read type. Please use " << "\'R\'" << " for sequencing reads or "
					"\'B\'" << " for barcode reads." << std::endl;
			return false;
		}

		temp.push_back(SequenceElement(temp.size(), (type == 'R') ? ++mates : 0, length));
		lenSum += length;

	}
  globalAlignmentSettings.set_seqs(temp);
  globalAlignmentSettings.set_mates(mates);

	if ( lenSum!=globalAlignmentSettings.get_cycles() ) {
		std::cerr << "Sum of defined reads does not equal the given number of cycles." << std::endl;
		return false;
	}

	return true;
}


bool parseBarcodeArgument(std::vector< std::string > barcodeArg ) {

	std::vector<uint16_t> barcode_lengths;

	for ( uint16_t seq_num = 0; seq_num < globalAlignmentSettings.get_seqs().size(); seq_num++ ) {

		// We are only interesed in Barcode sequences
		if ( !globalAlignmentSettings.getSeqById(seq_num).isBarcode() )
			continue;

		barcode_lengths.push_back( globalAlignmentSettings.getSeqById(seq_num).length );
	}

  std::vector<std::vector<std::string> > barcodeVector;
	for ( auto barcode = barcodeArg.begin(); barcode != barcodeArg.end(); ++barcode) {

		std::string valid_chars = seq_chars + "-";
		for(CountType i = 0; i != (*barcode).length(); i++){
			char c = (*barcode)[i];
			if ( valid_chars.find(c) == std::string::npos )
				return false;
		}

		std::vector<std::string> fragments;
		split(*barcode, '-', fragments);

		// check validity of barcode
		if ( barcode_lengths.size() != fragments.size())
			return false;

		for ( uint16_t num = 0; num != fragments.size(); num++ ) {
			if ( fragments[num].length() != barcode_lengths[num] ) {
				return false;
			}
		}

		// push back the fragments vector
		barcodeVector.push_back(fragments);
	}
	globalAlignmentSettings.set_barcodeVector(barcodeVector);

	return true;
}
