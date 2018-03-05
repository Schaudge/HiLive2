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
	po::options_description general("GENERAL OPTIONS");
	general.add_options()
	    		("help,h", "Print this help message and exit")
				("license,l", "Print licensing information and exit");
	return general;
}

po::options_description BuildIndexArgumentParser::positional_options() {
	po::options_description parameters("REQUIRED OPTIONS");
	parameters.add_options()
	    		("input,i", po::value<std::string>()->required(), "Reference genome(s) in (multi-)FASTA format. [REQUIRED]")
	    		("out-prefix,o", po::value<std::string>()->required(), "Output file prefix. Several files with the same prefix will be created. [REQUIRED]")
				;
	return parameters;
}

po::options_description BuildIndexArgumentParser::build_options() {
	po::options_description options("OTHER OPTIONS");
	options.add_options()
				("do-not-convert-spaces", po::bool_switch(&do_not_convert_spaces)->default_value(false), "Do not convert all spaces in reference ids to underscores [Default: converting is on]")
				("trim-after-space", po::bool_switch(&trim_ids)->default_value(false), "Trim all reference ids after first space [Default: false]")
				;
	return options;
}

void BuildIndexArgumentParser::init_help(po::options_description visible_options) {

	std::stringstream help_message;

	help_message << "Copyright (c) 2015-2018, Martin S. Lindner and the HiLive contributors. See CONTRIBUTORS for more info." << std::endl;
	help_message << "All rights reserved" << std::endl << std::endl;
	help_message << "HiLive is open-source software. Check with --license for details." << std::endl << std::endl;

	help_message << "Usage: " << std::endl << "  hilive-build [options]" << std::endl << std::endl;
	help_message << "Example command: " << std::endl << "  hilive-build --input hg19.fa --out-prefix ./index/hg19" << std::endl << std::endl;

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

	// File prefix of the index
	index_name = vm["outfile"].as<std::string>();

	return true;
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

	po::variables_map vm;
	try {
		// parse arguments
		po::store(po::command_line_parser(argc, argv).
				options(cmdline_options).run(), vm);
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

	  return 0;
}

//|||||||||||||||||||||||||||||||||||||||||||||||||||||||
//-----HiLiveArgumentParser------------------------------
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||

po::options_description HiLiveArgumentParser::general_options() {
	po::options_description general("GENERAL OPTIONS");
	general.add_options()
	    		  ("help,h", "Print this help message and exit.\n")
				  ("license,l", "Print license information and exit.\n")
				  ("config,c", po::value<std::string>(), "Path to a config file. Config file is in .ini format. Duplicate keys are not permitted. Instead, use comma-separated lists. Parameters obtained from the command line are prioritized over settings made in the config file.\n\nExample for a config.ini:\n   bcl-dir=./BaseCalls\n   lanes=1\n   out-cycle=50,100\n")
				  ("runinfo", po::value<std::string>(), "Path to Illumina's runInfo.xml file. If specified, read lengths, lane count and tile count are automatically set in accordance with the sequencing run. Parameters obtained from the command line or config file are prioritized over settings obtained from the runInfo.xml.\n")
				  ("continue", po::value<CountType>(), "Continue an interrupted HiLive run from a specified cycle. We strongly recommend to load the config file that was automatically created for the original run to continue with identical settings. This config file (hilive_config.ini) can be found in the temporary directory specified with --temp-dir.\n")
				  ;
	return general;
}

po::options_description HiLiveArgumentParser::sequencing_options() {
	po::options_description sequencing("SEQUENCING OPTIONS");
	sequencing.add_options()
					("bcl-dir,b", po::value<std::string>(), "Illumina's BaseCalls directory which contains the sequence information of the reads.\n")
					("lanes,l", po::value<std::string>(), "Specify the lanes to be considered for read alignment. [Default: 1-8]\n")
					("tiles,t", po::value<std::string>(), "Specify the tiles to be considered for read alignment. [Default: [1-2][1-3][01-16] (96 tiles)]\n")
					("reads,r", po::value<std::string>(), "Length and types of the read segments. Each segment is either a read ('R') or a barcode ('B'). Please give the segments in the correct order as they are produced by the sequencing machine. [REQUIRED]\n\nExample:\n   --reads 101R,8B,8B,101R\nspecifies paired-end sequencing with 2x101bp reads and 2x8bp barcodes.\n")
					("barcodes,B", po::value<std::string>(), "Barcode(s) of the sample(s) to be considered for read alignment. Barcodes must match the barcode length(s) as specified with --reads. Delimit different segments of the same barcodes by '-' and different barcodes by ','. [Default: All barcodes]\n\nExample:\n   -b ACCG-ATTG,ATGT-TGAC\nfor two different barcodes of length 2x4bp.\n")
					;
	return sequencing;
}

po::options_description HiLiveArgumentParser::report_options() {
	po::options_description report("REPORT OPTIONS");
	report.add_options()
					("out-dir,o", po::value<std::string>(), "Path to the directory that is used for the output files. The directory will be created if it does not exist. [Default: ./out]\n")
					("out-format,f", po::value<std::string>(), "Format of the output files. Currently, SAM and BAM format are supported. [Default: BAM]\n")
					("out-cycles,O", po::value<std::string>(), "Cycles for that alignment output is written. The respective temporary files are kept. [Default: write only after the last cycle]\n")
					("out-mode,M", po::value<std::string>(), "The output mode. [Default: ANYBEST]\n[ALL|A]: Report all found alignments.\n[BESTN#|N#]: Report the # best found alignments.\n[ALLBEST|H]: Report all found alignments with the best score.\n[ANYBEST|B]: Report one best alignment.\n[UNIQUE|U]: Report only unique alignments.\n")
					("report-unmapped", po::bool_switch(), "Activate reporting unmapped reads. [Default: false]\n")
					("extended-cigar", po::bool_switch(), "Activate extended CIGAR format for the alignment output files ('=' for matches and 'X' for mismatches instead of using 'M' for both). [Default: false]\n")
					("force-resort", po::bool_switch(), "Always sort temporary alignment files before writing output. Existing sorted align files are overwritten. This is only necessary if the temp directory is used more than once for new alignments. In general, this is not recommended for most applications. [Default: false (only sort if no sorted files exist)]\n")
					("max-softclip-ratio", po::value<float>(), "Maximal relative length of the front softclip (only relevant during output) [Default: 0.2]\n\nFurther explanation:\nHiLive uses an approach that requires one exact match of a k-mer at the beginning of an alignment. This can lead to unaligned regions at the beginning of the read which we report as 'softclips'. With this parameter, you can control the maximal length of this region.\n")
					;
	return report;
}

po::options_description HiLiveArgumentParser::alignment_options() {
	po::options_description alignment("ALIGNMENT OPTIONS");
	alignment.add_options()
					("index,i", po::value<std::string>(), "Path to the HiLive index. Please use the executable 'hilive-build' to create a new HiLive index that is delivered with this program. The index consists of several files with the same prefix. Please include the file prefix when specifying the index location.\n")
					("align-mode,m", po::value<std::string>(), "Alignment mode to balance speed and accuracy [fast|balanced|accurate]. This selected mode automatically sets other parameters. Individually configured parameters are prioritized over settings made by selecting an alignment mode. [Default: balanced]\n")
					("anchor-length", po::value<CountType>(), "Length of the alignment anchor (or initial seed) [Default: set by the selected alignment mode]\n")
					("error-interval", po::value<CountType>(), "The interval to tolerate more errors during alignment (low=accurate; great=fast). [Default: 'anchor-length'/2]\n")
					("seeding-interval", po::value<CountType>(), "The interval to create new seeds (low=accurate; great=fast). [Default: 'anchor-length'/2]\n")
					("barcode-errors", po::value<std::string>(), "The number of errors that are tolerated for the barcode segments. A single value can be provided to be applied for all barcode segments. Alternatively, the value can be set for each segment individually. [Default: 1]\n\nExample:\n   --barcode-errors 2 [2 errors for all barcode segments]\n   --barcode-errors 2,1 [2 errors for the first, 1 error for the second segment]\n")
					("align-undetermined-barcodes", po::bool_switch()->default_value(false), "Align all barcodes. Reads with a barcode that don't match one of the barcodes specified with '--barcodes' will be reported as undetermined. [Default: false]\n")
					("min-basecall-quality", po::value<CountType>(), "Minimum basecall quality for a nucleotide to be considered as a match [Default: 1 (everything but N-calls)]\n")
					("keep-invalid-sequences", po::bool_switch(), "Keep sequences of invalid reads, i.e. with unconsidered barcode or filtered by the sequencer. This option must be activated to report unmapped reads. [Default: false]\n")
					;
	return alignment;
}

po::options_description HiLiveArgumentParser::scoring_options() {
	po::options_description scoring("SCORING OPTIONS");
	scoring.add_options()
					("min-as,s", po::value<ScoreType>(), "Minimum alignment score. [Default: Set automatically based on the alignment mode and match/mismatch scores]\n")
					("match-score", po::value<CountType>(), "Score for a match. [Default: 0]\n")
					("mismatch-penalty", po::value<CountType>(), "Penalty for a mismatch. [Default: 6]\n")
					("insertion-opening-penalty", po::value<CountType>(), "Penalty for insertion opening. [Default: 5]\n")
					("insertion-extension-penalty", po::value<CountType>(), "Penalty for insertion extension. [Default: 3]\n")
					("deletion-opening-penalty", po::value<CountType>(), "Penalty for deletion opening. [Default: 5]\n")
					("deletion-extension-penalty", po::value<CountType>(), "Penalty for deletion extension. [Default: 3]\n")
					("max-gap-length", po::value<CountType>(), "Maximal permitted consecutive gap length. Increasing this parameter may lead to highly increased runtime! [Default: 3]\n")
					("softclip-opening-penalty", po::value<float>(), "Penalty for softclip opening (only relevant during output). [Default: 'mismatch-penalty']\n")
					("softclip-extension-penalty", po::value<float>(), "Penalty for softclip extension (only relevant during output). [Default: 'mismatch-penalty'/'anchor-length']\n")
					;
	return scoring;

}

po::options_description HiLiveArgumentParser::technical_options() {
	 po::options_description technical("TECHNICAL OPTIONS");
	 technical.add_options()
	    	        ("temp-dir", po::value<std::string>(), "Temporary directory to store the alignment files and hilive_config.ini. [Default: ./temp]\n")
					("keep-files,k", po::value<std::string>(), "Keep intermediate alignment files for these cycles. The last cycle is always kept. [Default: Keep files of output cycles]\n\nFurther Explanations:\nHiLive comes with a separated executable 'hilive-out'. This executable can be used to produce alignment files in SAM or BAM format from existing temporary files. Thus, output can only be created for cycles for that keeping the temporary alignment files is activated. Temporary alignemnt files are also needed if an interrupted run is continued with the '--continue' parameter.\n")
					("keep-all-files,K", po::bool_switch(), "Keep all intermediate alignment files. This option may lead to huge disk space requirements. [Default: false]\n")
	        		("block-size", po::value<std::string>(), "Block size for the alignment input/output stream in Bytes. Append 'K' or 'M' to specify in Kilobytes or Megabytes, respectively. [Default: 64M]\n\nExample:\n   --block-size 1024 [1024 bytes]\n   --block-size 64K [64 Kilobytes]\n   --block-size 64M [64 Megabytes]\n")
					("compression", po::value<uint16_t>(), "Compression of temporary alignment files. [Default: LZ4]\n0: no compression.\n1: Deflate (smaller).\n2: LZ4 (faster).\n")
					("num-threads,n", po::value<CountType>(), "Number of threads to spawn (including output threads). [Default: 1]\n")
					("num-out-threads,N", po::value<CountType>(), "Maximum number of threads to use for output. More threads may be used for output automatically if threads are idle. [Default: 'num-threads'/2]\n")
					;
	 return technical;
}

void HiLiveArgumentParser::init_help(po::options_description visible_options) {

	std::stringstream help_message;

	help_message << "Copyright (c) 2015-2018, Martin S. Lindner and the " << std::endl << "HiLive contributors. See CONTRIBUTORS for more info." << std::endl;
	help_message << "All rights reserved" << std::endl << std::endl;
	help_message << "HiLive is open-source software. Check with --license for details." << std::endl << std::endl;

	help_message << "Usage: " << std::endl << "  hilive [options]" << std::endl << std::endl;
	help_message << "Example command: " << std::endl << "  hilive --bcl-dir ./BaseCalls --index ./reference/hg19 --reads 101R" << std::endl << std::endl;
	help_message << "REQUIRED OPTIONS:" << std::endl;
	help_message << "  --bcl-dir             Illumina's BaseCalls directory which contains the sequence information of the reads." << std::endl;
	help_message << "  --index               Path to the HiLive index." << std::endl;
	help_message << "  --reads               Length and types of the read segments." << std::endl << std::endl;;
	help_message << "Required options might be specified either on the command line or in the config file." << std::endl;

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

	// TODO: REACTIVATE!

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
//				if ( sequences.size() > 0 )
//					runInfo_settings.add_child("settings.sequences", getXMLnode_vector(sequences));
				runInfo_settings.put("settings.cycles", num_cycles);

	            if (ptree_Run.count("FlowcellLayout")!=0) {

	            	ptree ptree_FlowcellLayout = ptree_Run.get_child("FlowcellLayout");

	            	// Get the lanes
	            	std::vector<uint16_t> lanes_vec(ptree_FlowcellLayout.get<unsigned>("<xmlattr>.LaneCount"));
	            	std::iota(lanes_vec.begin(), lanes_vec.end(), 1);
//	            	runInfo_settings.add_child("settings.lanes", getXMLnode_vector(lanes_vec));

	            	// Get the tiles
	            	std::vector<uint16_t> tiles_vec = flowcell_layout_to_tile_numbers(
	            			ptree_FlowcellLayout.get<unsigned>("<xmlattr>.SurfaceCount"),
							ptree_FlowcellLayout.get<unsigned>("<xmlattr>.SwathCount"),
							ptree_FlowcellLayout.get<unsigned>("<xmlattr>.TileCount") );
//	            	runInfo_settings.add_child("settings.tiles", getXMLnode_vector(tiles_vec));
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
		set_option<std::string>("reads", "settings.sequences", "", &AlignmentSettings::set_read_structure);

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
		set_option<std::string>("out-cycles", "settings.out.cycles", "" + globalAlignmentSettings.get_cycles(), &AlignmentSettings::set_output_cycles);
		set_option<bool>("extended-cigar", "settings.out.extended_cigar", false, &AlignmentSettings::set_extended_cigar);

		set_option<ScoreType>("min-as", "settings.out.min_as", ScoreType(getMaxPossibleScore(globalAlignmentSettings.getSeqByMate(1).length) - 3*globalAlignmentSettings.get_mismatch_penalty()), &AlignmentSettings::set_min_as); // TODO: change default

		if ( cmd_settings.at("keep-all-files").as<bool>() ) {
			std::vector<CountType>keep_all_files (globalAlignmentSettings.get_cycles());
			std::iota(keep_all_files.begin(), keep_all_files.end(), 1);
			globalAlignmentSettings.set_keep_aln_files(keep_all_files);
		} else {
			set_option<std::string>("keep-files", "settings.technical.keep_aln_files", "", &AlignmentSettings::set_keep_aln_files);
		}
		set_option<std::string>("lanes", "settings.lanes", to_string(all_lanes()), &AlignmentSettings::set_lanes);
		set_option<std::string>("tiles", "settings.tiles", to_string(all_tiles()), &AlignmentSettings::set_tiles);
		set_option<CountType>("seeding-interval", "settings.align.seeding_interval", globalAlignmentSettings.get_anchor_length()/2, &AlignmentSettings::set_seeding_interval);


		set_option<bool>("force-resort", "settings.out.force-resort", false, &AlignmentSettings::set_force_resort);

		set_option<std::string>("out-mode", "settings.mode", "ANYBEST", &AlignmentSettings::set_mode);
		set_option<CountType>("min-basecall-quality", "settings.align.min_qual", 1, &AlignmentSettings::set_min_qual);

		std::vector<std::string> barcode_sequences_default;
		set_option<std::vector<std::string>>("barcodes", "settings.barcodes.sequences", barcode_sequences_default, &AlignmentSettings::set_barcodes);

		set_option<std::string>("barcode-errors", "settings.barcodes.errors", "2", &AlignmentSettings::set_barcode_errors);

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

	help_message << "Usage: " << std::endl << "  hilive-out --config /path/to/settings/file [options]" << std::endl << std::endl;
	help_message << "REQUIRED OPTIONS:" << std::endl;
	help_message << "  --config              Path to a HiLive config file (in general, this should be 'hilive_config.ini' which is created in the temp directory of the respective run)" << std::endl;
	help_message << std::endl << "All parameters can be set as for the HiLive main program." << std::endl;
	help_message << "By default, only output files for the last cycle are produced." << std::endl;
	help_message << "Use the --out-cycles parameter to declare different cycle numbers (will only work if --keep-files or --out-cycles was activated for the cycle when running HiLive)" << std::endl;

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
