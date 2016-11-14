#include "argument_parser.h"
namespace po = boost::program_options;


int parseCommandLineArguments(AlignmentSettings & settings, std::string license, int argc, char const ** argv)
{
    po::options_description general("General");
    general.add_options()
        ("help,h", "Print this help message and exit")
        ("license", "Print licensing information and exit");

    po::options_description parameters("Parameters");
    parameters.add_options()
        ("BC_DIR", po::value<std::string>(&settings.root)->required(), "Illumina BaseCalls directory")
        ("INDEX", po::value<std::string>(&settings.index_fname)->required(), "Path to k-mer index")
        ("CYCLES", po::value<CountType>(&settings.cycles)->required(), "Number of cycles")
        ("OUTDIR", po::value<std::string>(&settings.out_dir), "Directory to store sam files in [Default: temporary or BaseCalls directory");

    po::options_description io_settings("IO settings");
    io_settings.add_options()
        ("temp", po::value<std::string>(&settings.temp_dir)->default_value(""), "Temporary directory for the alignment files [Default: use BaseCalls directory]")
        //("bam,B", po::bool_switch(&settings.write_bam)->default_value(false), "Create BAM files instead of SAM files [Default: false]")
        ("keep-files,k", po::bool_switch(&settings.keep_aln_files)->default_value(false), "Keep intermediate alignment files [Default: false]")
        ("lanes,l", po::value< std::vector<uint16_t> >()->multitoken()->composing(), "Select lane [Default: all lanes]")
        ("tiles,t", po::value< std::vector<uint16_t> >()->multitoken()->composing(), "Select tile numbers [Default: all tiles]")
        ("barcodes,b", po::value< std::vector<std::string> >()->multitoken()->composing(), "Enumerate barcodes (must have same length) for demultiplexing, e.g. -b AGGATC -b CCCTTT [Default: no demultiplexing]")
    	("barcode-errors,E", po::value< std::vector<uint16_t> >()->multitoken()->composing(), "Enumerate the number of tolerated errors (only SNPs) for each barcode fragment, e.g. -E 2 2 [Default: 1 per fragment]")
		("reads,r", po::value< std::vector<std::string> >()->multitoken()->composing(), "Enumerate read lengths and type. Example: -r 101R 8B 8B 101R equals paired-end sequencing with 2x101bp reads and 2x8bp barcodes.");

    po::options_description alignment("Alignment settings");
    alignment.add_options()
        ("min-errors,e", po::value<CountType>(&settings.min_errors)->default_value(2), "Number of errors tolerated in read alignment [Default: 2]")
        ("all-best-hit,H", po::bool_switch()->default_value(false), "Report all of the best alignments for each read")
        ("any-best-hit", po::bool_switch(), "Report one of the best alignments for each read (default)")
        ("all-best-n-scores,N", po::value<CountType>(&settings.best_n), "Report all alignments of the N best alignment scores for each read")
        ("all-hits,A", po::bool_switch()->default_value(false), "Report all valid alignments for each read")
        ("disable-ohw-filter", po::bool_switch(&settings.discard_ohw)->default_value(true), "disable the One-Hit Wonder filter [Default: false]")
        ("start-ohw", po::value<CountType>(&settings.start_ohw)->default_value(K_HiLive+5), "First cycle to apply One-Hit Wonder filter [Default: K+5]")
        ("window,w", po::value<DiffType>(&settings.window)->default_value(5), "Set the window size to search for alignment extension, i.e. maximum total insertion/deletion size [Default: 5]")
        ("min-quality", po::value<CountType>(&settings.min_qual)->default_value(1), "Minimum allowed basecall quality [Default: 1]");

    po::options_description technical("Technical settings");
    technical.add_options()
        ("block-size", po::value<uint64_t>()->default_value(64*1024*1024), "Block size for the alignment input/output stream in Bytes. Use -K or -M to specify in Kilobytes or Megabytes")
        (",K", po::bool_switch()->default_value(false), "Interpret the block-size argument as Kilobytes instead of Bytes")
        (",M", po::bool_switch()->default_value(false), "Interpret the block-size argument as Megabytes instead of Bytes")
        ("compression,c", po::value<uint8_t>(&settings.compression_format)->default_value(2), "Compress alignment files. 0: no compression (default) 1: Deflate (smaller) 2: LZ4 (faster)")
        ("num-threads,n", po::value<int>(), "Number of threads to spawn [Default: all available]");

    po::options_description cmdline_options;
    cmdline_options.add(general).add(parameters).add(io_settings).add(alignment).add(technical);

    po::options_description visible_options;
    visible_options.add(general).add(io_settings).add(alignment).add(technical);

    std::stringstream help_message;
    help_message << "HiLive v"<< HiLive_VERSION_MAJOR << "." << HiLive_VERSION_MINOR << " - Realtime Alignment of Illumina Reads" << std::endl;
    help_message << "Copyright (c) 2015, Martin S. Lindner" << std::endl;
    help_message << "HiLive is open-source software. Check with --license for details." << std::endl << std::endl;
    help_message << "Fixed k-mer size: " << K_HiLive << std::endl << std::endl;
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

    if (vm.count("OUTDIR"))
        settings.out_dir = vm["OUTDIR"].as<std::string>();
    else {
        if (settings.temp_dir == "") 
            settings.out_dir = settings.root;
        else
            settings.out_dir = settings.temp_dir;
    }

    // Parse read lengths and types. If read argument is missing, init for single-end and non-barcoded.
    if ( vm.count("reads") && !parseReadsArgument(settings, vm["reads"].as< std::vector<std::string> >()) )
    	return -1;
    else if ( !vm.count("reads") ) {
    	settings.seqs.push_back(SequenceElement(0,1,settings.cycles));
    	settings.mates = 1;
    }

    if (vm.count("lanes"))
        settings.lanes = vm["lanes"].as< std::vector<uint16_t> >();
    else
        settings.lanes = all_lanes();

    if (vm.count("tiles"))
        settings.tiles = vm["tiles"].as< std::vector<uint16_t> >();
    else
        settings.tiles = all_tiles();



    settings.barcodeVector.clear();
    if (vm.count("barcodes")) {
        settings.barcodeVector = vm["barcodes"].as< std::vector<std::string> >();
        if( !parseBarcodeArgument(settings, vm["barcodes"].as< std::vector<std::string> >()) ) {
        	std::cerr << "Parsing error: Invalid barcode(s) detected. Please ensure that you used \'-\' "
        			"as duplex delimiter and that all barcodes have the correct length. Only use A,C,G and T as bases!" << std::endl;
        	return -1;
        }
    }

    if ( settings.multiBarcodeVector.size() != 0 ) {
    	if ( vm.count("barcode-errors") ) {
    		settings.barcode_errors = vm["barcode-errors"].as< std::vector<uint16_t> >();
    		if ( settings.multiBarcodeVector[0].size() != settings.barcode_errors.size() ) {
    			std::cerr << "Parsing error: Number of barcode errors does not equal the number of barcodes." << std::endl;
    			return -1;
    		}
    	} else {
        	for ( uint16_t i = 0; i < settings.multiBarcodeVector[0].size(); i++ ) {
        		settings.barcode_errors.push_back(1);
        	}
    	}
    }

    if (vm["all-hits"].as<bool>()) {
        // all hits: disable other modes
        settings.any_best_hit_mode = false;
        settings.all_best_hit_mode = false;
        settings.all_best_n_scores_mode = false;
    } else if (vm["all-best-hit"].as<bool>()) {
        // enable all-best-hit mode and disable others
        settings.any_best_hit_mode = false;
        settings.all_best_hit_mode = true;
        settings.all_best_n_scores_mode = false;
    } else if (vm.count("all-best-n-scores")) {
        // enable any-best-n mode and get parameter
        settings.any_best_hit_mode = false;
        settings.all_best_hit_mode = false;
        settings.all_best_n_scores_mode = true;
        settings.best_n = vm["all-best-n-scores"].as<CountType>();
    } else { // the default behaviour
        // enable any-best-hit mode and disable others
        settings.any_best_hit_mode = true;
        settings.all_best_hit_mode = false;
        settings.all_best_n_scores_mode = false;
    }
        
    if (vm["-M"].as<bool>())
        settings.block_size = vm["block-size"].as<uint64_t>()*1024*1024;
    else if (vm["-K"].as<bool>())
        settings.block_size = vm["block-size"].as<uint64_t>()*1024;
    else
        settings.block_size = vm["block-size"].as<uint64_t>();

    if (vm.count("num-threads")) 
        settings.num_threads = vm["num-threads"].as<int>();
    else { // the default case, meaning as much as physically useful
        uint32_t n_cpu = std::thread::hardware_concurrency();
        if (n_cpu > 1)
            settings.num_threads = n_cpu;
        else
            settings.num_threads = 1;
    }


    // check paths and file names
    if (!file_exists(settings.index_fname)){
        std::cerr << "Input error: Could not find k-mer index file " << settings.index_fname << std::endl;
        return -1;
    }

    std::size_t found = settings.root.find("BaseCalls");
    if (!(found != std::string::npos && found >= settings.root.size()-10)) {
        std::cerr << "Warning: BaseCalls directory seems to be invalid: " << settings.root << std::endl;
    } 

    if (!is_directory(settings.root)){
        std::cerr << "Input error: Could not find BaseCalls directory " << settings.root << std::endl;
        return -1;
    }

    for ( uint16_t ln : settings.lanes ) {
        std::string ln_dir = settings.root;
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


    // Report the basic settings
    std::cout << "Running HiLive with       " << settings.num_threads << " thread(s)." << std::endl;
    std::cout << "BaseCalls directory:      " << settings.root << std::endl;
    if (settings.temp_dir != "") {
        std::cout << "Temporary directory:      " << settings.temp_dir << std::endl;
    }
    //if (!settings.write_bam)
    std::cout << "SAM output directory:     " << settings.out_dir << std::endl;
    //else
        //std::cout << "BAM output directory:     " << settings.out_dir << std::endl;
    std::cout << "Lanes:                    ";
    for ( uint16_t ln : settings.lanes )
        std::cout << ln << " ";
    std::cout << std::endl;
    std::cout << "K-mer index:              " << settings.index_fname << std::endl;
    std::cout << "Read lengths:             ";
    std::string barcode_suffix;
    for ( uint16_t read = 0; read != settings.seqs.size(); read ++) {
    	std::cout << settings.getSeqById(read).length;
    	barcode_suffix = settings.getSeqById(read).isBarcode() ? "B" : "R";
    	std::cout << barcode_suffix << " ";
    }
    std::cout << std::endl;
    std::cout << "Mapping error:            " << settings.min_errors << std::endl;
    if (settings.any_best_hit_mode) 
        std::cout << "Mapping mode:             Any-Best-Hit-Mode" << std::endl;
    else if (settings.all_best_hit_mode) 
        std::cout << "Mapping mode:             All-Best-Hit-Mode" << std::endl;
    else if (settings.all_best_n_scores_mode) 
        std::cout << "Mapping mode:             All-Best-N-Scores-Mode with N=" << settings.best_n << std::endl;
    else
        std::cout << "Mapping mode:             All-Hits-Mode" << std::endl;
    std::cout << std::endl;

    return 0;
}

bool parseReadsArgument(AlignmentSettings & settings, std::vector< std::string > readsArg){

	CountType lenSum = 0;
	CountType length = 0;
	std::string length_string = "";
	char type;

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

		settings.seqs.push_back(SequenceElement(settings.seqs.size(), (type == 'R') ? ++settings.mates : 0, length));
		lenSum += length;

	}

	if ( lenSum!=settings.cycles ) {
		std::cerr << "Sum of defined reads does not equal the given number of cycles." << std::endl;
		return false;
	}

	return true;
}

bool parseBarcodeArgument(AlignmentSettings & settings, std::vector< std::string > barcodeArg ) {

	std::vector<uint16_t> barcode_lengths;

	for ( uint16_t seq_num = 0; seq_num < settings.seqs.size(); seq_num++ ) {

		// We are only interesed in Barcode sequences
		if ( !settings.getSeqById(seq_num).isBarcode() )
			continue;

		barcode_lengths.push_back( settings.getSeqById(seq_num).length );
	}

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
		settings.multiBarcodeVector.push_back(fragments);

	}

	return true;

}
