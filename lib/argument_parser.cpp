#include "argument_parser.h"

namespace po = boost::program_options;

//|||||||||||||||||||||||||||||||||||||||||||||||||||||||
//-----ArgumentParser------------------------------------
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||

ArgumentParser::ArgumentParser(int argC, char const ** argV):argc(argC),argv(argV){}

//|||||||||||||||||||||||||||||||||||||||||||||||||||||||
//-----BuildIndexArgumentParser--------------------------
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||

po::options_description BuildIndexArgumentParser::general_options() {
	po::options_description general("General");
	general.add_options()
	    		("help,h", "Print this help message and exit")
				("license", "Print licensing information and exit");
	return general;
}

po::options_description BuildIndexArgumentParser::positional_options() {
	po::options_description parameters("Required");
	parameters.add_options()
	    		("INPUT", po::value<std::string>()->required(), "Reference genomes in (multi-) FASTA format.");
	return parameters;
}

po::options_description BuildIndexArgumentParser::build_options() {
	po::options_description options("Options");
	options.add_options()
	    		("outfile,o", po::value<std::string>(), "Set output file name [Default: INPUT.kix]")
				("do-not-convert-spaces", po::bool_switch(&do_not_convert_spaces)->default_value(false), "Do not convert all spaces in reference ids to underscores [Default: converting is on]")
				("trim-after-space", po::bool_switch(&trim_ids)->default_value(false), "Trim all reference ids after first space [Default: false]");
	return options;
}

void BuildIndexArgumentParser::init_help(po::options_description visible_options) {

	std::stringstream help_message;

	help_message << "Copyright (c) 2015-2017, Martin S. Lindner and the HiLive contributors. See CONTRIBUTORS for more info." << std::endl;
	help_message << "All rights reserved" << std::endl << std::endl;
	help_message << "HiLive is open-source software. Check with --license for details." << std::endl << std::endl;
	help_message << "Usage: " << std::endl << "  hilive-build INPUT [options]" << std::endl << std::endl;
	help_message << "Required: " << std::endl;
	help_message << "  INPUT                 Reference genomes in (multi-) FASTA format." << std::endl;

	help_message << visible_options;

	help = help_message.str();
}

bool BuildIndexArgumentParser::set_positional_variables(po::variables_map vm) {

	// Name of the input fasta file
	fasta_name = vm["INPUT"].as<std::string>();

	// Check if input file exists
	if ( !file_exists(fasta_name) ){
		std::cerr << "Input error: Could not find input file " << fasta_name << std::endl;
		return false;
	}

	return true;
}

bool BuildIndexArgumentParser::set_build_variables(po::variables_map vm) {

	// Name of the ouput .kix file
	if (vm.count("outfile")) {
		index_name = vm["outfile"].as<std::string>();
	} else {
		index_name = fasta_name + std::string(".fmix");
	}

	return true;
}

void BuildIndexArgumentParser::report() {

//	std::cout << "K-mer weight:        " << (uint16_t) globalAlignmentSettings.get_kmer_weight() << std::endl;
//	std::cout << "K-mer span:          " << (uint16_t) globalAlignmentSettings.get_kmer_span() << std::endl;
//	std::cout << "K-mer gap positions: ";
//
//	if ( globalAlignmentSettings.get_kmer_gaps().size() > 0 ) {
//		for ( auto pos : globalAlignmentSettings.get_kmer_gaps() ) {
//			if ( pos != *(globalAlignmentSettings.get_kmer_gaps().begin()) )
//				std::cout << ",";
//			std::cout << (uint16_t) pos;
//		}
//		std::cout << std::endl;
//	} else {
//		std::cout << "-" << std::endl;
//	}
//	std::cout << std::endl;

}

int BuildIndexArgumentParser::parseCommandLineArguments() {

	po::options_description gen_opt = general_options();
	po::options_description pos_opt = positional_options();
	po::options_description build_opt = build_options();

	po::options_description cmdline_options;
	cmdline_options.add(pos_opt).add(gen_opt).add(build_opt);

	po::options_description visible_options;
	visible_options.add(gen_opt).add(build_opt);

	init_help(visible_options);

	po::positional_options_description p;
	p.add("INPUT", 1);

	po::variables_map vm;
	try {
		// parse arguments
		po::store(po::command_line_parser(argc, argv).
				options(cmdline_options).positional(p).run(), vm);
		// first check if -h or --help was called
		if (vm.count("help")) {
			printHelp();
			return 1;
	    }
	    // first check if --license was called
	    if (vm.count("license")) {
	      printLicense();
	      return 1;
	    }

	    // then check arguments
	    po::notify(vm);
	  }
	  catch ( po::required_option& e ) {
	    std::cerr << "Missing Parameter: " << e.what() << std::endl;
	    return -1;
	  }
	  catch( po::error& e) {
	    std::cerr << "Error while parsing command line options: " << e.what() << std::endl;
	    return -1;
	  }

	  if ( !set_positional_variables(vm) ) {
		  return -1;
	  }

	  if ( !set_build_variables(vm) ) {
		  return -1;
	  }

	  report();

	  return 0;
}

//|||||||||||||||||||||||||||||||||||||||||||||||||||||||
//-----HiLiveArgumentParser------------------------------
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||

po::options_description HiLiveArgumentParser::general_options() {
	po::options_description general("--- GENERAL OPTIONS ---");
	general.add_options()
	    		  ("help", "Print this help message and exit")
				  ("license", "Print licensing information and exit")
				  ("config", po::value<std::string>(), "Load settings from file. If command line arguments are given additionally, they are prefered.")
				  ("continue", po::value<CountType>(), "Continue an interrupted HiLive run from a specified cycle. We strongly recommend to load the settings from the previous run using the -s option.")
				  ;
	return general;
}

po::options_description HiLiveArgumentParser::sequencing_options() {
	po::options_description sequencing("--- SEQUENCING OPTIONS ---");
	sequencing.add_options()
					("bcl-dir", po::value<std::string>(), "Illumina BaseCalls directory")
					("runinfo", po::value<std::string>(), "Path to runInfo.xml for parsing read and index lengths [Default (if activated): BC_DIR/../../RunInfo.xml]")
					("lanes", po::value< std::vector<uint16_t> >()->multitoken()->composing(), "Select lane [Default: all lanes]")
					("tiles", po::value< std::vector<uint16_t> >()->multitoken()->composing(), "Select tile numbers [Default: all tiles]")
					("reads", po::value< std::vector<std::string> >()->multitoken()->composing(), "Enumerate read lengths and type. Example: -r 101R 8B 8B 101R equals paired-end sequencing with 2x101bp reads and 2x8bp barcodes. Overwrites information of runInfo.xml. [Default: single end reads without barcodes]")
					("barcodes", po::value< std::vector<std::string> >()->multitoken()->composing(), "Enumerate barcode(s) of your sample(s), e.g. -b AGGATC -b CCCTTT [Default: Align all barcodes and report them as 'undetermined']")
					;
	return sequencing;
}

po::options_description HiLiveArgumentParser::report_options() {
	po::options_description report("--- REPORT OPTIONS ---");
	report.add_options()
					("out-dir", po::value<std::string>(), "Directory for output files [Default: ./out]")
					("out-format", po::value<std::string>(), "Format of the output files. SAM and BAM format are currently supported [Default: BAM]")
					("out-cycles", po::value<std::vector<CountType>>()->multitoken()->composing(), "Cycles for alignment output. The respective temporary files are kept. [Default: last cycle]")
					("out-mode", po::value<std::string>(), "Output mode. [ALL|A]: Report all alignments; [BESTN#|N#]: Report alignments of the best # scores; "
							"[ALLBEST|H]: Report all alignments with the best score (similar to N1); [ANYBEST|B]: Report one best alignment (default)")
					("report-unmapped", po::bool_switch(), "Report unmapped reads [Default: false]")
					("extended-cigar", po::bool_switch(), "Activate extended CIGAR format (= and X instead of only M) in output files [Default: false]")
					("force-resort", po::bool_switch(), "If set, temporary files are always sorted before output. Existing sorted align files are overwritten [Default: false]")
					;
	return report;
}

po::options_description HiLiveArgumentParser::alignment_options() {
	po::options_description alignment("--- ALIGNMENT OPTIONS ---");
	alignment.add_options()
					("index", po::value<std::string>(), "Path to the HiLive index.")
					("align-mode", po::value<std::string>(), "Alignment mode to balance speed and accuracy automatically. [fast|balanced|accurate] [Default: balanced]")
					("anchor-length", po::value<CountType>(), "Set the anchor length manually [Default: 12]")
					("error-interval", po::value<CountType>(), "Set the interval to allow more errors (low=accurate; great=fast) [Default: anchor-length/2]")
					("seeding-interval", po::value<CountType>(), "Set the interval to create new seeds (low=accurate; great=fast) [Default: error-interval]")
					("barcode-errors", po::value< std::vector<uint16_t> >()->multitoken()->composing(), "Enumerate the number of tolerated errors (only SNPs) for each barcode fragment, e.g. -E 2 2 [Default: 1 per fragment]")
					("align-undetermined-barcodes", po::bool_switch()->default_value(false), "Align and output all barcodes [Default: false]")
					("min-basecall-quality", po::value<CountType>(), "Minimum allowed basecall quality [Default: 1]")
					("keep-invalid-sequences", po::bool_switch(), "Keep sequences of invalid reads, i.e. with unconsidered barcode or filtered by the sequencer (needed for reporting those unmapped reads) [Default: false]")
					;
	return alignment;
}

po::options_description HiLiveArgumentParser::scoring_options() {
	po::options_description scoring("--- SCORING OPTIONS ---");
	scoring.add_options()
					("min-as", po::value<ScoreType>(), "Minimum alignment score [Default: Depends on score model]")
					("match-score", po::value<CountType>(), "Score for a match [Default: 0]")
					("mismatch-penalty", po::value<CountType>(), "Penalty for a mismatch [Default: 6]")
					("insertion-opening-penalty", po::value<CountType>(), "Penalty for insertion opening [Default: 5]")
					("insertion-extension-penalty", po::value<CountType>(), "Penalty for insertion extension [Default: 3]")
					("deletion-opening-penalty", po::value<CountType>(), "Penalty for deletion opening [Default: 5]")
					("deletion-extension-penalty", po::value<CountType>(), "Penalty for deletion extension [Default: 3]")
					("max-gap-length", po::value<CountType>(), "Maximal gap length. High influence on runtime depending on the scoring scheme [Default: 3]")
					("softclip-opening-penalty", po::value<float>(), "Penalty for softclip opening (only relevant during output!) [Default: mismatch-penalty]")
					("softclip-extension-penalty", po::value<float>(), "Penalty for softclip extension (only relevant during output!) [Default: mismatch-penalty/error-rate]")
					("max-softclip-ratio", po::value<float>(), "Maximal relative length of the front softclip (only relevant during output!) [Default: 0.2]")
					;
	return scoring;

}

po::options_description HiLiveArgumentParser::technical_options() {
	 po::options_description technical("--- TECHNICAL OPTIONS ---");
	 technical.add_options()
	    	        ("temp-dir", po::value<std::string>(), "Temporary directory for the alignment files [Default: ./temp]")
					("keep-files,k", po::value<std::vector<CountType>>()->multitoken()->composing(), "Keep intermediate alignment files for these cycles. The last cycle is always kept. [Default: None]")
					("keep-all-files,K", po::bool_switch(), "Keep all intermediate alignment files [Default: false]")
	        		("block-size", po::value<std::string>(), "Block size for the alignment input/output stream in Bytes. Append 'K' or 'M' to specify in Kilobytes or Megabytes, respectively (e.g. '--block-size 64M' for 64 Megabytes)")
					("compression,c", po::value<uint16_t>(), "Compress alignment files. 0: no compression 1: Deflate (smaller) 2: LZ4 (faster; default)")
					("num-threads,n", po::value<CountType>(), "Number of threads to spawn [Default: all available]")
					("num-out-threads,N", po::value<CountType>(), "Maximum number of threads to use for output if threads are not idle [Default: half of -n]")
					;
	 return technical;
}

void HiLiveArgumentParser::init_help(po::options_description visible_options) {

	std::stringstream help_message;

	help_message << "Copyright (c) 2015-2017, Martin S. Lindner and the HiLive contributors. See CONTRIBUTORS for more info." << std::endl;
	help_message << "All rights reserved" << std::endl << std::endl;
	help_message << "HiLive is open-source software. Check with --license for details." << std::endl << std::endl;

	help_message << "Usage: " << std::endl << "  hilive [options]" << std::endl << std::endl;
	help_message << "REQUIRED OPTIONS:" << std::endl;
	help_message << "  --bcl-dir                Illumina BaseCalls directory of the sequencing run to analyze" << std::endl;
	help_message << "  --index                  Path to the HiLive index" << std::endl;
	help_message << "  --reads                  Sequencing cycles specified for each mate" << std::endl;

	help_message << visible_options;

	help = help_message.str();
}

bool HiLiveArgumentParser::checkPaths() {
//	if (!file_exists(globalAlignmentSettings.get_index_fname())){
//		std::cerr << "Input error: Could not find k-mer index file " << globalAlignmentSettings.get_index_fname() << std::endl;
//		return false;
//	}

	std::size_t found = globalAlignmentSettings.get_root().find("BaseCalls");
	if (!(found != std::string::npos && found >= globalAlignmentSettings.get_root().size()-10)) {
		std::cerr << "Warning: BaseCalls directory seems to be invalid: " << globalAlignmentSettings.get_root() << std::endl;
	}

	if (!is_directory(globalAlignmentSettings.get_root())){
		std::cerr << "Input error: Could not find BaseCalls directory " << globalAlignmentSettings.get_root() << std::endl;
		return false;
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
			return false;
		}
	}

	if ( !is_directory(globalAlignmentSettings.get_temp_dir())) {
		boost::filesystem::create_directories(globalAlignmentSettings.get_temp_dir());
	}

	if ( !is_directory(globalAlignmentSettings.get_out_dir())) {
		boost::filesystem::create_directories(globalAlignmentSettings.get_out_dir());
	}

	return true;
}

void HiLiveArgumentParser::report() {
    std::cout << "Running HiLive with       " << globalAlignmentSettings.get_num_threads() << " thread(s)." << std::endl;
    std::cout << "BaseCalls directory:      " << globalAlignmentSettings.get_root() << std::endl;
    if (globalAlignmentSettings.get_temp_dir() != "") {
        std::cout << "Temporary directory:      " << globalAlignmentSettings.get_temp_dir() << std::endl;
    }
    if ( globalAlignmentSettings.get_output_format() == OutputFormat::SAM )
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
    std::cout << "Min. alignment score:     " << globalAlignmentSettings.get_min_as() << std::endl;
    if (globalAlignmentSettings.get_any_best_hit_mode())
        std::cout << "Mapping mode:             Any-Best-Hit-Mode" << std::endl;
    else if (globalAlignmentSettings.get_all_best_hit_mode())
        std::cout << "Mapping mode:             All-Best-Hit-Mode" << std::endl;
    else if (globalAlignmentSettings.get_all_best_n_scores_mode())
        std::cout << "Mapping mode:             All-Best-N-Scores-Mode with N=" << globalAlignmentSettings.get_best_n() << std::endl;
    else if (globalAlignmentSettings.get_unique_hit_mode())
        std::cout << "Mapping mode:             Unique-Hits-Mode" << std::endl;
    else
        std::cout << "Mapping mode:             All-Hits-Mode" << std::endl;

    std::cout << "Anchor length:            " << globalAlignmentSettings.get_anchor_length() << std::endl;

	if ( globalAlignmentSettings.get_start_cycle() > 1 ) {
		std::cout << std::endl;
		std::cout << "----- CONTINUE RUN FROM CYCLE " << cmd_settings.at("continue").as<CountType>() << " -----" << std::endl;
	}

    std::cout << std::endl;
}

bool HiLiveArgumentParser::parseRunInfo(po::variables_map vm) {

	if ( ! vm.count("runinfo"))
		return false;

	boost::property_tree::ptree tree;
	read_xml(tree, vm["runinfo"].as<std::string>());

	using boost::property_tree::ptree;

	if (!tree.empty() && tree.count("RunInfo")!=0) {
		ptree ptree_RunInfo = tree.get_child("RunInfo");

		if (ptree_RunInfo.count("Run")!=0) {
			ptree ptree_Run = ptree_RunInfo.get_child("Run");

			if (ptree_Run.count("Reads")!=0) {

				ptree ptree_Reads = ptree_Run.get_child("Reads");

				// Get the sequence structure and total number of cycles
				std::vector<std::string> sequences;
				CountType num_cycles = 0;
				for (const auto &read : ptree_Reads) {
					std::string sequence = "";
					sequence += read.second.get<std::string>("<xmlattr>.NumCycles");
					sequence += read.second.get<std::string>("<xmlattr>.IsIndexedRead") == "N" ? "R" : "B";
					sequences.push_back(sequence);
					num_cycles += read.second.get<unsigned>("<xmlattr>.NumCycles");
				}
				if ( sequences.size() > 0 )
					runInfo_settings.add_child("settings.sequences", getXMLnode_vector(sequences));
				runInfo_settings.put("settings.cycles", num_cycles);

	            if (ptree_Run.count("FlowcellLayout")!=0) {

	            	ptree ptree_FlowcellLayout = ptree_Run.get_child("FlowcellLayout");

	            	// Get the lanes
	            	std::vector<uint16_t> lanes_vec(ptree_FlowcellLayout.get<unsigned>("<xmlattr>.LaneCount"));
	            	std::iota(lanes_vec.begin(), lanes_vec.end(), 1);
	            	runInfo_settings.add_child("settings.lanes", getXMLnode_vector(lanes_vec));

	            	// Get the tiles
	            	std::vector<uint16_t> tiles_vec = flowcell_layout_to_tile_numbers(
	            			ptree_FlowcellLayout.get<unsigned>("<xmlattr>.SurfaceCount"),
							ptree_FlowcellLayout.get<unsigned>("<xmlattr>.SwathCount"),
							ptree_FlowcellLayout.get<unsigned>("<xmlattr>.TileCount") );
	            	runInfo_settings.add_child("settings.tiles", getXMLnode_vector(tiles_vec));
	            }
			}
		}
	}
	return true;
}

int HiLiveArgumentParser::parseCommandLineArguments() {

	this->set_required_parameters();

	// Init general options
	po::options_description gen_opt = general_options();
	po::variables_map vm;

	// Init all other options
	po::options_description sequence_opt = sequencing_options();
	po::options_description report_opt = report_options();
	po::options_description align_opt = alignment_options();
	po::options_description score_opt = scoring_options();
	po::options_description tech_opt = technical_options();

	// All command line options
    po::options_description cmdline_options;
    cmdline_options.add(gen_opt).add(sequence_opt).add(report_opt).add(align_opt).add(score_opt).add(tech_opt);

    // Options visible in the help
    po::options_description visible_options;
    visible_options.add(gen_opt).add(sequence_opt).add(report_opt).add(align_opt).add(score_opt).add(tech_opt);

    init_help(visible_options);

    // First parameter iteration for general options (includes help, license and input settings file.
    try {
        po::store(po::command_line_parser(argc, argv).options(gen_opt).allow_unregistered().run(), vm);

        // first check if -h or --help was called
        if (vm.count("help")) {
        	printHelp();
        	return 1;
        }

        // first check if --license was called
        if (vm.count("license")) {
        	printLicense();
        	return 1;
        }

        vm.notify();

        // Load input settings if exist
        if ( vm.count("config")) {
        	if ( file_exists(vm["config"].as<std::string>() ) ) {
        		std::ifstream config_file(vm["config"].as<std::string>());
        		po::store(po::parse_config_file(config_file, cmdline_options, false), input_settings);
        	} else {
        		throw file_not_exist_error(vm["config"].as<std::string>());
        	}
        } else if ( isRequired("config") ) {
        	throw po::required_option("config");
        }

        if ( vm.count("runinfo") ) {
        	if ( ! parseRunInfo(vm) ) {
                std::cerr << "Error while parsing Run Info file: " << vm["runinfo"].as<std::string>() << std::endl;
                return -1;

        	}
        } else if ( isRequired("runinfo") ) {
        	throw po::required_option("runinfo");
        }
    }
    catch( po::error& e) {
           std::cerr << "Error while parsing command line options: " << std::endl << e.what() << std::endl;
           return -1;
    }

    po::positional_options_description p;
    p.add("bcl-dir", 1);
    p.add("index", 1);
    p.add("reads", 1);

    // Parse all command line arguments to cmd_settings
    try {
        // parse arguments
        po::store(po::parse_command_line(argc, argv, cmdline_options), cmd_settings);

        // then check arguments
        po::notify(cmd_settings);
        po::notify(input_settings);
    }
    catch ( po::required_option& e ) {
        std::cerr << "Missing Parameter: " << e.what() << std::endl;
        return -1;  
    }
    catch( po::error& e) {
        std::cerr << "Error while parsing command line options: " << e.what() << std::endl;
        return -1;  
    } 

    // Set all options from command line and input files
    if ( !set_options() ) {
    	return -1;
    }

    if ( !checkPaths() ) {
    	return -1;
    }

    // Report the basic settings
    report();

    return 0;
}

bool HiLiveArgumentParser::set_options() {

	try {

		// Set continue cycle if given by the user
		if ( cmd_settings.count("continue") ) {
			globalAlignmentSettings.set_start_cycle(cmd_settings.at("continue").as<CountType>());
		} else {
			globalAlignmentSettings.set_start_cycle(1);
		}

		// Set positional arguments
		set_option<std::string>("bcl-dir", "settings.paths.root", "", &AlignmentSettings::set_root);

		set_option<std::string>("index", "settings.paths.index", "", &AlignmentSettings::set_index_fname);

		set_option<std::string>("out-dir", "settings.paths.out_dir", "./out", &AlignmentSettings::set_out_dir);


		CountType default_anchor_length = set_mode();

		// Set arguments that are required for setting other arguments
		set_option<CountType>("anchor-length", "settings.align.anchor", default_anchor_length, &AlignmentSettings::set_anchor_length);
		set_option<CountType>("error-interval", "settings.align.error_interval", globalAlignmentSettings.get_anchor_length()/2, &AlignmentSettings::set_error_rate);


		// Set I/O options
		set_option<std::string>("temp-dir", "settings.paths.temp_dir", "./temp", &AlignmentSettings::set_temp_dir);
		set_option<std::string>("out-format", "settings.out.bam", "BAM", &AlignmentSettings::set_output_format);
		set_option<bool>("keep-invalid-sequences", "settings.out.keep_all_sequences", false, &AlignmentSettings::set_keep_all_sequences);
		set_option<bool>("report-unmapped", "settings.out.report_unmapped", false, &AlignmentSettings::set_report_unmapped);

		// Set read structure
		set_option<std::vector<std::string>>("reads", "settings.sequences", std::vector<std::string>(), &AlignmentSettings::set_read_structure);

		// Scoring scheme
		set_option<CountType>("match-score", "settings.scores.match_score", 0, &AlignmentSettings::set_match_score);
		set_option<CountType>("mismatch-penalty", "settings.scores.mismatch_penalty", 6, &AlignmentSettings::set_mismatch_penalty);
		set_option<CountType>("insertion-opening-penalty", "settings.scores.insertion_opening_penalty", 5, &AlignmentSettings::set_insertion_opening_penalty);
		set_option<CountType>("deletion-opening-penalty", "settings.scores.deletion_opening_penalty", 5, &AlignmentSettings::set_deletion_opening_penalty);
		set_option<CountType>("insertion-extension-penalty", "settings.scores.insertion_extension_penalty", 3, &AlignmentSettings::set_insertion_extension_penalty);
		set_option<CountType>("deletion-extension-penalty", "settings.scores.deletion_extension_penalty", 3, &AlignmentSettings::set_deletion_extension_penalty);
		set_option<CountType>("max-gap-length", "settings.scores.max_gap_length", 3, &AlignmentSettings::set_max_gap_length);
		set_option<float>("softclip-opening-penalty", "settings.scores.softclip_opening_penalty", float(globalAlignmentSettings.get_mismatch_penalty()), &AlignmentSettings::set_softclip_opening_penalty);
		set_option<float>("softclip-extension-penalty", "settings.scores.softclip_extension_penalty", float(globalAlignmentSettings.get_mismatch_penalty()) / globalAlignmentSettings.get_anchor_length(), &AlignmentSettings::set_softclip_extension_penalty);
		set_option<float>("max-softclip-ratio", "settings.scores.max_softclip_ratio", .2f, &AlignmentSettings::set_max_softclip_ratio);


		// Alignment options
		std::vector<CountType> output_cycles = {globalAlignmentSettings.get_cycles()};
		set_option<std::vector<CountType>>("output-cycles", "settings.out.cycles", output_cycles, &AlignmentSettings::set_output_cycles);
		set_option<bool>("extended-cigar", "settings.out.extended_cigar", false, &AlignmentSettings::set_extended_cigar);

		set_option<ScoreType>("min-as", "settings.out.min_as", ScoreType(getMaxPossibleScore(globalAlignmentSettings.getSeqByMate(1).length) - 3*globalAlignmentSettings.get_mismatch_penalty()), &AlignmentSettings::set_min_as); // TODO: change default

		if ( cmd_settings.at("keep-all-files").as<bool>() ) {
			std::vector<CountType>keep_all_files (globalAlignmentSettings.get_cycles());
			std::iota(keep_all_files.begin(), keep_all_files.end(), 1);
			globalAlignmentSettings.set_keep_aln_files(keep_all_files);
		} else {
			set_option<std::vector<CountType>>("keep-files", "settings.technical.keep_aln_files", std::vector<CountType>(), &AlignmentSettings::set_keep_aln_files);
		}
		set_option<std::vector<uint16_t>>("lanes", "settings.lanes", all_lanes(), &AlignmentSettings::set_lanes);
		set_option<std::vector<uint16_t>>("tiles", "settings.tiles", all_tiles(), &AlignmentSettings::set_tiles);
		set_option<CountType>("seeding-interval", "settings.align.seeding_interval", globalAlignmentSettings.get_anchor_length()/2, &AlignmentSettings::set_seeding_interval);


		set_option<bool>("force-resort", "settings.out.force-resort", false, &AlignmentSettings::set_force_resort);

		set_option<std::string>("out-mode", "settings.mode", "ANYBEST", &AlignmentSettings::set_mode);
		set_option<CountType>("min-basecall-quality", "settings.align.min_qual", 1, &AlignmentSettings::set_min_qual);

		std::vector<std::string> barcode_sequences_default;
		set_option<std::vector<std::string>>("barcodes", "settings.barcodes.sequences", barcode_sequences_default, &AlignmentSettings::set_barcodes);

		std::vector<CountType> barcode_errors_default = {2};
		set_option<std::vector<CountType>>("barcode-errors", "settings.barcodes.errors", barcode_errors_default, &AlignmentSettings::set_barcode_errors);

		set_option<bool>("align-undetermined-barcodes", "settings.barcodes.keep_all", false, &AlignmentSettings::set_keep_all_barcodes);

		// Set technical options
		set_option<std::string>("block-size", "settings.technical.block_size", "64M", &AlignmentSettings::set_block_size);
		set_option<uint16_t>("compression", "settings.technical.compression_format", 2, &AlignmentSettings::set_compression_format);

		CountType n_cpu = std::thread::hardware_concurrency();
		CountType n_threads_default = 1;
		if (n_cpu > 1)
			n_threads_default = std::min( n_cpu, CountType( globalAlignmentSettings.get_lanes().size() * globalAlignmentSettings.get_tiles().size() ) ) ;
		set_option<CountType>("num-threads", "settings.technical.num_threads", n_threads_default, &AlignmentSettings::set_num_threads);
		set_option<CountType>("num-out-threads", "settings.technical.num_out_threads", globalAlignmentSettings.get_num_threads()/2, &AlignmentSettings::set_num_out_threads);

	} catch ( std::exception & ex ) {
		std::cerr << "Error while parsing options: " << std::endl << ex.what() << std::endl;
		return false;
	}
	return true;
}

CountType HiLiveArgumentParser::set_mode() {

	CountType default_anchor_length = 15;
	uint64_t genome_size = 0;

	KixRun* tempIdx = new KixRun();
	tempIdx->load_seqlengths(globalAlignmentSettings.get_index_fname());
	for (CountType i=0; i<tempIdx->getNumSequences(); i++ ) {
		genome_size += 2*tempIdx->getSeqLengths()[i];
	}

	// relative number of reads matching a reference of given length randomly (in theory, not in biology)
	float expectation_value = .0025f;

	CountType balanced_anchor_length = ( std::log(float(genome_size) / expectation_value ) / std::log(4) );

	// Only return default value if no mode is set.
	if ( !cmd_settings.count("alignment-mode") )
		return balanced_anchor_length;

	char mode = std::toupper(cmd_settings["alignment-mode"].as<std::string>()[0]);

	// Accurate
	if( mode=='A' ) {
		default_anchor_length = std::floor(0.833f * balanced_anchor_length);

	// Balanced
	} else if ( mode == 'B' ) {
		default_anchor_length = balanced_anchor_length;

	// Fast
	} else if ( mode == 'F' ) {
		default_anchor_length = std::ceil(1.166f * balanced_anchor_length);
	}

	else {
		throw po::invalid_option_value ("--alignment-mode " + cmd_settings["alignment-mode"].as<std::string>());
	}

	// Insert default values to variables map if not already set.
	if ( !cmd_settings.count("anchor-length") )
		cmd_settings.insert(std::make_pair("anchor-length", po::variable_value(default_anchor_length, true)));

	if ( !cmd_settings.count("error-interval") )
		cmd_settings.insert(std::make_pair("error-interval", po::variable_value(CountType(default_anchor_length/2), true)));

	if ( !cmd_settings.count("seeding-interval") )
		cmd_settings.insert(std::make_pair("seeding-interval", po::variable_value(CountType(default_anchor_length/2), true)));


	return default_anchor_length;
}

//|||||||||||||||||||||||||||||||||||||||||||||||||||||||
//-----HiLiveOutArgumentParser--------------------------
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||

void HiLiveOutArgumentParser::init_help(po::options_description visible_options) {

	std::stringstream help_message;

	help_message << "Copyright (c) 2015-2017, Martin S. Lindner and the HiLive contributors. See CONTRIBUTORS for more info." << std::endl;
	help_message << "All rights reserved" << std::endl << std::endl;
	help_message << "HiLive is open-source software. Check with --license for details." << std::endl << std::endl;

	help_message << "Usage: " << std::endl << "  hilive-out --settings /path/to/settings/file [options]" << std::endl << std::endl;
	help_message << "Required:" << std::endl;
	help_message << "  settings              Path to a HiLive settings file (by default, the file is in the temp directory of the respective run)" << std::endl;
	help_message << std::endl << "All parameters can be set as for the HiLive main program." << std::endl;
	help_message << "By default, only output files for the last cycle are produced." << std::endl;
	help_message << "Use the --output-cycles parameter to declare different cycle numbers (will only work if --keep-files was activated for the respective HiLive run)" << std::endl;

	help_message << visible_options;

	help = help_message.str();
}

void HiLiveOutArgumentParser::report() {
	if (globalAlignmentSettings.get_temp_dir() != "") {
	        std::cout << "Temporary directory:      " << globalAlignmentSettings.get_temp_dir() << std::endl;
	}
	if ( globalAlignmentSettings.get_output_format() == OutputFormat::SAM )
		std::cout << "SAM output directory:     " << globalAlignmentSettings.get_out_dir() << std::endl;
	else
		std::cout << "BAM output directory:     " << globalAlignmentSettings.get_out_dir() << std::endl;
	std::cout << "Lanes:                    ";
	for ( uint16_t ln : globalAlignmentSettings.get_lanes() )
		std::cout << ln << " ";
	std::cout << std::endl;
	std::cout << "K-mer index:              " << globalAlignmentSettings.get_index_fname() << std::endl;
	std::cout << "Total Read lengths:       ";
	std::string barcode_suffix;
	for ( uint16_t read = 0; read != globalAlignmentSettings.get_seqs().size(); read ++) {
		std::cout << globalAlignmentSettings.getSeqById(read).length;
		barcode_suffix = globalAlignmentSettings.getSeqById(read).isBarcode() ? "B" : "R";
		std::cout << barcode_suffix << " ";
	}
	std::cout << std::endl;
	std::cout << "Min. Alignment Score:     " << globalAlignmentSettings.get_min_as() << std::endl;
	if (globalAlignmentSettings.get_any_best_hit_mode())
		std::cout << "Mapping mode:             Any-Best-Hit-Mode" << std::endl;
	else if (globalAlignmentSettings.get_all_best_hit_mode())
		std::cout << "Mapping mode:             All-Best-Hit-Mode" << std::endl;
	else if (globalAlignmentSettings.get_all_best_n_scores_mode())
		std::cout << "Mapping mode:             All-Best-N-Scores-Mode with N=" << globalAlignmentSettings.get_best_n() << std::endl;
    else if (globalAlignmentSettings.get_unique_hit_mode())
        std::cout << "Mapping mode:             Unique-Hits-Mode" << std::endl;
	else
		std::cout << "Mapping mode:             All-Hits-Mode" << std::endl;
	std::cout << "Output Cycles:            ";
	for ( auto cycle : globalAlignmentSettings.get_output_cycles() ) {
		std::cout << cycle << " ";
	}
	std::cout << std::endl;
	std::cout << std::endl;
}
