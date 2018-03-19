#include "argument_parser.h"

namespace po = boost::program_options;

//|||||||||||||||||||||||||||||||||||||||||||||||||||||||
//-----ArgumentParser------------------------------------
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||

ArgumentParser::ArgumentParser(int argC, char const ** argV):argc(argC),argv(argV){}

std::string ArgumentParser::select_prioritized_parameter( std::vector<std::string> parameters ) {

	// Check command line for parameters
	for ( auto & param : parameters ) {
		if ( cmd_settings.count ( param ) ) {
			return param;
		}
	}

	// Check config file for parameters
	for ( auto & param : parameters ) {
		if ( config_file_settings.count ( param ) ) {
			return param;
		}
	}

	// Check runInfo for parameters
	for ( auto & param : parameters ) {
		if ( config_file_settings.count ( param ) ) {
			return param;
		}
	}

	return parameters.front();
}


//|||||||||||||||||||||||||||||||||||||||||||||||||||||||
//-----BuildIndexArgumentParser--------------------------
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||

po::options_description BuildIndexArgumentParser::general_options() {
	po::options_description general("GENERAL OPTIONS", default_line_width, default_line_width*.9f);
	general.add_options()
	    		("help,h", "Print this help message and exit")
				("license,l", "Print licensing information and exit");
	return general;
}

po::options_description BuildIndexArgumentParser::positional_options() {
	po::options_description parameters("REQUIRED OPTIONS", default_line_width, default_line_width*.9f);
	parameters.add_options()
	    		("input,i", po::value<std::string>()->required(), "Reference genome(s) in (multi-)FASTA format. [REQUIRED]")
	    		("out-prefix,o", po::value<std::string>()->required(), "Output file prefix. Several files with the same prefix will be created. [REQUIRED]")
				;
	return parameters;
}

po::options_description BuildIndexArgumentParser::build_options() {
	po::options_description options("OTHER OPTIONS", default_line_width, default_line_width*.9f);
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

	help_message << visible_options;

	help = help_message.str();
}

bool BuildIndexArgumentParser::set_positional_variables(po::variables_map vm) {

	// Name of the input fasta file
	fasta_name = vm["input"].as<std::string>();

	// Check if input file exists
	if ( !file_exists(fasta_name) ){
		std::cerr << "Input error: Could not find input file " << fasta_name << std::endl;
		return false;
	}

	// File prefix of the index
	index_name = vm["out-prefix"].as<std::string>();

	return true;
}

int BuildIndexArgumentParser::parseCommandLineArguments() {

	po::options_description gen_opt = general_options();
	po::options_description pos_opt = positional_options();
	po::options_description build_opt = build_options();

	po::options_description cmdline_options;
	cmdline_options.add(pos_opt).add(gen_opt).add(build_opt);

	po::options_description visible_options;
	visible_options.add(gen_opt).add(pos_opt).add(build_opt);

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
	po::options_description general("GENERAL OPTIONS", default_line_width, default_line_width*.9f);
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
	po::options_description sequencing("SEQUENCING OPTIONS", default_line_width, default_line_width*.9f);
	sequencing.add_options()
					("bcl-dir,b", po::value<std::string>(), "Illumina's BaseCalls directory which contains the sequence information of the reads.\n")
					("lanes,l", po::value<std::string>(), "Specify the lanes to be considered for read alignment. [Default: 1-8]\n")
					("tiles,t", po::value<std::string>(), "Specify the tiles to be considered for read alignment. [Default: [1-2][1-3][01-16] (96 tiles)]\n")
					("max-tile,T", po::value<CountType>(), "Specify the highest tile number. The tile numbers will be computed by this number, considering the correct surface count, swath count and tile count for Illumina sequencing.\nThis parameter serves as a shortcut for --tiles.\n\nExample:\n   --max-tile 2216\nwill activate all tiles in [1-2][1-2][01-16].\n")
					("reads,r", po::value<std::string>(), "Length and types of the read segments. Each segment is either a read ('R') or a barcode ('B'). Please give the segments in the correct order as they are produced by the sequencing machine. [REQUIRED]\n\nExample:\n   --reads 101R,8B,8B,101R\nspecifies paired-end sequencing with 2x101bp reads and 2x8bp barcodes.\n")
					("barcodes,B", po::value<std::string>(), "Barcode(s) of the sample(s) to be considered for read alignment. Barcodes must match the barcode length(s) as specified with --reads. Delimit different segments of the same barcodes by '-' and different barcodes by ','. [Default: All barcodes]\n\nExample:\n   -b ACCG-ATTG,ATGT-TGAC\nfor two different barcodes of length 2x4bp.\n")
					("run-id", po::value<std::string>(), "ID of the sequencing run. Should be obtained from runInfo.xml.\n")
					("flowcell-id", po::value<std::string>(), "ID of the flowcell. Should be obtained from runInfo.xml.\n")
					("instrument-id", po::value<std::string>(), "ID of the sequencing machine. Should be obtained from runInfo.xml.\n")
					;
	return sequencing;
}

po::options_description HiLiveArgumentParser::report_options() {
	po::options_description report("REPORT OPTIONS", default_line_width, default_line_width*.9f);
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
	po::options_description alignment("ALIGNMENT OPTIONS", default_line_width, default_line_width*.9f);
	alignment.add_options()
					("index,i", po::value<std::string>(), "Path to the HiLive index. Please use the executable 'hilive-build' to create a new HiLive index that is delivered with this program. The index consists of several files with the same prefix. Please include the file prefix when specifying the index location.\n")
					("align-mode,m", po::value<std::string>(), "Alignment mode to balance speed and accuracy [fast|balanced|accurate]. This selected mode automatically sets other parameters. Individually configured parameters are prioritized over settings made by selecting an alignment mode. [Default: balanced]\n")
					("anchor-length", po::value<CountType>(), "Length of the alignment anchor (or initial seed) [Default: set by the selected alignment mode]\n")
					("error-interval", po::value<CountType>(), "The interval to tolerate more errors during alignment (low=accurate; great=fast). [Default: 'anchor-length'/2]\n")
					("seeding-interval", po::value<CountType>(), "The interval to create new seeds (low=accurate; great=fast). [Default: 'anchor-length'/2]\n")
					("barcode-errors", po::value<std::string>(), "The number of errors that are tolerated for the barcode segments. A single value can be provided to be applied for all barcode segments. Alternatively, the value can be set for each segment individually. [Default: 1]\n\nExample:\n   --barcode-errors 2 [2 errors for all barcode segments]\n   --barcode-errors 2,1 [2 errors for the first, 1 error for the second segment]\n")
					("align-undetermined-barcodes", po::bool_switch(), "Align all barcodes. Reads with a barcode that don't match one of the barcodes specified with '--barcodes' will be reported as undetermined. [Default: false]\n")
					("min-basecall-quality", po::value<CountType>(), "Minimum basecall quality for a nucleotide to be considered as a match [Default: 1 (everything but N-calls)]\n")
					("keep-invalid-sequences", po::bool_switch(), "Keep sequences of invalid reads, i.e. with unconsidered barcode or filtered by the sequencer. This option must be activated to report unmapped reads. [Default: false]\n")
					;
	return alignment;
}

po::options_description HiLiveArgumentParser::scoring_options() {
	po::options_description scoring("SCORING OPTIONS", default_line_width, default_line_width*.9f);
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
	 po::options_description technical("TECHNICAL OPTIONS", default_line_width, default_line_width*.9f);
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
	help_message << "  -b [--bcl-dir]        Illumina's BaseCalls directory which contains the sequence information of the reads." << std::endl;
	help_message << "  -i [--index]          Path to the HiLive index." << std::endl;
	help_message << "  -r [--reads]          Length and types of the read segments." << std::endl << std::endl;;
	help_message << "Required options might be specified either on the command line or in the config file." << std::endl;

	help_message << visible_options;

	help = help_message.str();
}

bool HiLiveArgumentParser::checkPaths() {

	std::size_t found = globalAlignmentSettings.get_root().find("BaseCalls");
	if (!(found != std::string::npos && found >= globalAlignmentSettings.get_root().size()-9)) {
		std::cerr << "Warning: BaseCalls directory seems to be invalid: " << globalAlignmentSettings.get_root() << std::endl;
	}

	if (!is_directory(globalAlignmentSettings.get_root())){
		std::cerr << "Input error: Could not find BaseCalls directory " << globalAlignmentSettings.get_root() << std::endl;
		return false;
	}

	for ( uint16_t ln : globalAlignmentSettings.get_lanes() ) {
		std::string ln_dir = globalAlignmentSettings.get_root() + "/L" + to_N_digits(ln, 3);
		if (!is_directory(ln_dir)){
			std::cout << "Waiting for directory of lane L" << to_N_digits(ln, 3) << " (" << ln_dir << ")..." << std::endl;
			while ( !is_directory(ln_dir) )
		        std::this_thread::sleep_for(std::chrono::seconds(60));
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
    	std::cout << globalAlignmentSettings.get_seq_by_id(read).length;
    	barcode_suffix = globalAlignmentSettings.get_seq_by_id(read).isBarcode() ? "B" : "R";
    	std::cout << barcode_suffix << " ";
    }
    std::cout << std::endl;
    std::cout << "Min. alignment score:     " << globalAlignmentSettings.get_min_as() << std::endl;
    std::cout << "Mapping mode:             " << to_string(globalAlignmentSettings.get_mode(), globalAlignmentSettings.get_best_n()) << std::endl;
    std::cout << "Anchor length:            " << globalAlignmentSettings.get_anchor_length() << std::endl;

	if ( globalAlignmentSettings.get_start_cycle() > 1 ) {
		std::cout << std::endl;
		std::cout << "----- CONTINUE RUN FROM CYCLE " << cmd_settings.at("continue").as<CountType>() << " -----" << std::endl;
	}

    std::cout << std::endl;
}

bool HiLiveArgumentParser::parseRunInfo(po::variables_map vm) {

	// Check for runinfo commandline parameter
	if ( ! vm.count("runinfo"))
		return false;

	using boost::property_tree::ptree;

	// Load the file
	ptree tree;
	read_xml(tree, vm["runinfo"].as<std::string>());

	// Try to obtain the run ID
	if ( tree.get_child_optional("RunInfo.Run.<xmlattr>.Id")) {
			runInfo_settings.insert(std::make_pair("run-id", boost::program_options::variable_value(tree.get<std::string>("RunInfo.Run.<xmlattr>.Id"), false)));
	}

	// Try to obtain the flowcell ID
	if ( tree.get_child_optional("RunInfo.Run.Flowcell") ) {
		runInfo_settings.insert(std::make_pair("flowcell-id", boost::program_options::variable_value(tree.get<std::string>("RunInfo.Run.Flowcell"), false)));
	}

	// Try to obtain the instrument ID
	if ( tree.get_child_optional("RunInfo.Run.Instrument")) {
		runInfo_settings.insert(std::make_pair("instrument-id", boost::program_options::variable_value(tree.get<std::string>("RunInfo.Run.Instrument"), false)));
	}

	 // Try to obtain the read segments
	if ( tree.get_child_optional("RunInfo.Run.Reads") ) {

		std::vector<std::string> sequences;
		for (const auto &read : tree.get_child("RunInfo.Run.Reads") ) {

			if ( !read.second.get_child_optional("<xmlattr>.NumCycles")
					|| !read.second.get_child_optional("<xmlattr>.IsIndexedRead") ) {
				throw std::runtime_error("Parsing error: Read information in runInfo file is not valid.");
			}

			// Get the segments
			std::string sequence = "";
			sequence += read.second.get<std::string>("<xmlattr>.NumCycles");
			sequence += read.second.get<std::string>("<xmlattr>.IsIndexedRead") == "N" ? "R" : "B";
			sequences.push_back(sequence);
		}

		// Store the read segments
		runInfo_settings.insert(std::make_pair("reads", boost::program_options::variable_value(join(sequences), false)));
	}

	// Try to obtain the Flowcell layout
	if ( tree.get_child_optional("RunInfo.Run.FlowcellLayout") ) {

		auto tree_FlowcellLayout = tree.get_child("RunInfo.Run.FlowcellLayout");

		// Store the lane count
		if ( tree_FlowcellLayout.get_child_optional("<xmlattr>.LaneCount")) {
			std::vector<uint16_t> lanes_vec(tree_FlowcellLayout.get<unsigned>("<xmlattr>.LaneCount"));
			std::iota(lanes_vec.begin(), lanes_vec.end(), 1);
			runInfo_settings.insert(std::make_pair("lanes", boost::program_options::variable_value(join(lanes_vec), false)));
		}

		// Store the tile count
		if ( tree_FlowcellLayout.get_child_optional("<xmlattr>.SurfaceCount")
				&& tree_FlowcellLayout.get_child_optional("<xmlattr>.SwathCount")
				&& tree_FlowcellLayout.get_child_optional("<xmlattr>.TileCount") ) {

			std::vector<uint16_t> tiles_vec = flowcell_layout_to_tile_numbers(
					tree_FlowcellLayout.get<unsigned>("<xmlattr>.SurfaceCount"),
					tree_FlowcellLayout.get<unsigned>("<xmlattr>.SwathCount"),
					tree_FlowcellLayout.get<unsigned>("<xmlattr>.TileCount") );
			runInfo_settings.insert(std::make_pair("tiles", boost::program_options::variable_value(join(tiles_vec), false)));
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
        		po::store(po::parse_config_file(config_file, cmdline_options, false), config_file_settings);
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
        po::notify(config_file_settings);
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

		// GENERAL OPTIONS

		set_option<CountType>("continue", 1, &AlignmentSettings::set_start_cycle);


		// SEQUENCING OPTIONS

		set_option<std::string>("bcl-dir", "", &AlignmentSettings::set_root);

		set_option<std::string>("lanes", join(all_lanes()), &AlignmentSettings::set_lanes);

		if ( select_prioritized_parameter( {"tiles", "max-tile"} ) == "tiles" )
			set_option<std::string>("tiles", join(all_tiles()), &AlignmentSettings::set_tiles);
		else
			set_option<CountType>("max-tile", 2316, &AlignmentSettings::set_max_tile);

		set_option<std::string>("reads", "", &AlignmentSettings::set_read_structure);

		set_option<std::string>("barcodes", "", &AlignmentSettings::set_barcodes);

		set_option<std::string>("run-id", "", &AlignmentSettings::set_run_id);

		set_option<std::string>("flowcell-id", "", &AlignmentSettings::set_flowcell_id);

		set_option<std::string>("instrument-id", "", &AlignmentSettings::set_instrument_id);


		// REPORT OPTIONS

		set_option<std::string>("out-dir", "./out", &AlignmentSettings::set_out_dir);

		set_option<std::string>("out-format", "BAM", &AlignmentSettings::set_output_format);

		set_option<std::string>("out-cycles", std::to_string(globalAlignmentSettings.get_cycles()), &AlignmentSettings::set_output_cycles);

		set_option<std::string>("out-mode", "ANYBEST", &AlignmentSettings::set_mode);

		set_option<bool>("report-unmapped", false, &AlignmentSettings::set_report_unmapped);

		set_option<bool>("extended-cigar", false, &AlignmentSettings::set_extended_cigar);

		set_option<bool>("force-resort", false, &AlignmentSettings::set_force_resort);

		set_option<float>("max-softclip-ratio", .2f, &AlignmentSettings::set_max_softclip_ratio);


		// ALIGNMENT OPTIONS

		set_option<std::string>("index", "", &AlignmentSettings::set_index_fname);

		CountType default_anchor_length = set_mode();
		set_option<CountType>("anchor-length", default_anchor_length, &AlignmentSettings::set_anchor_length);

		set_option<CountType>("error-interval", globalAlignmentSettings.get_anchor_length()/2, &AlignmentSettings::set_error_rate);

		set_option<CountType>("seeding-interval", globalAlignmentSettings.get_anchor_length()/2, &AlignmentSettings::set_seeding_interval);

		set_option<std::string>("barcode-errors", "2", &AlignmentSettings::set_barcode_errors);

		set_option<bool>("align-undetermined-barcodes", false, &AlignmentSettings::set_keep_all_barcodes);

		set_option<CountType>("min-basecall-quality", 1, &AlignmentSettings::set_min_qual);

		set_option<bool>("keep-invalid-sequences", false, &AlignmentSettings::set_keep_all_sequences);


		// SCORING OPTIONS

		set_option<CountType>("match-score", 0, &AlignmentSettings::set_match_score);

		set_option<CountType>("mismatch-penalty", 6, &AlignmentSettings::set_mismatch_penalty);

		set_option<CountType>("insertion-opening-penalty", 5, &AlignmentSettings::set_insertion_opening_penalty);

		set_option<CountType>("deletion-opening-penalty", 5, &AlignmentSettings::set_deletion_opening_penalty);

		set_option<CountType>("insertion-extension-penalty", 3, &AlignmentSettings::set_insertion_extension_penalty);

		set_option<CountType>("deletion-extension-penalty", 3, &AlignmentSettings::set_deletion_extension_penalty);

		set_option<CountType>("max-gap-length", 3, &AlignmentSettings::set_max_gap_length);

		set_option<float>("softclip-opening-penalty", float(globalAlignmentSettings.get_mismatch_penalty()), &AlignmentSettings::set_softclip_opening_penalty);

		set_option<float>("softclip-extension-penalty", float(globalAlignmentSettings.get_mismatch_penalty()) / globalAlignmentSettings.get_anchor_length(), &AlignmentSettings::set_softclip_extension_penalty);

		// 3% error rate by default
		ScoreType min_as_default = getMaxPossibleScore(globalAlignmentSettings.get_seq_by_mate(1).length) - ( float(globalAlignmentSettings.get_seq_by_mate(1).length / 100.0f) * 3.0f * getMaxSingleErrorPenalty());
		set_option<ScoreType>("min-as", min_as_default, &AlignmentSettings::set_min_as);


		// TECHNICAL OPTIONS

		set_option<std::string>("temp-dir", "./temp", &AlignmentSettings::set_temp_dir);

		if ( select_prioritized_parameter( {"keep-all-files", "keep-files"} ) == "keep-all-files" ) {
			globalAlignmentSettings.set_keep_all_aln_files();
		} else {
			set_option<std::string>("keep-files", "", &AlignmentSettings::set_keep_aln_files);
		}

		set_option<std::string>("block-size", "64M", &AlignmentSettings::set_block_size);

		set_option<uint16_t>("compression", 2, &AlignmentSettings::set_compression_format);

		set_option<CountType>("num-threads", 1, &AlignmentSettings::set_num_threads);

		set_option<CountType>("num-out-threads", globalAlignmentSettings.get_num_threads()/2, &AlignmentSettings::set_num_out_threads);

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

	help_message << "Copyright (c) 2015-2018, Martin S. Lindner and the HiLive contributors. See CONTRIBUTORS for more info." << std::endl;
	help_message << "All rights reserved" << std::endl << std::endl;
	help_message << "HiLive is open-source software. Check with --license for details." << std::endl << std::endl;

	help_message << "Usage: " << std::endl << "  hilive-out --config /path/to/config/file [options]" << std::endl << std::endl;
	help_message << "REQUIRED OPTIONS:" << std::endl;
	help_message << "  -c [ --config ]" << std::endl << "        Path to a HiLive config file (in general, this should be" << std::endl << "        'hilive_config.ini' which is created in the temp directory of the" << std::endl << "        respective run)" << std::endl;
	help_message << std::endl << "All parameters can be set as for the HiLive main program." << std::endl;
	help_message << "By default, output cycles are the same as specified for the original HiLive run." << std::endl;
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
		std::cout << globalAlignmentSettings.get_seq_by_id(read).length;
		barcode_suffix = globalAlignmentSettings.get_seq_by_id(read).isBarcode() ? "B" : "R";
		std::cout << barcode_suffix << " ";
	}
	std::cout << std::endl;
	std::cout << "Min. Alignment Score:     " << globalAlignmentSettings.get_min_as() << std::endl;
	std::cout << "Mapping mode:             " << to_string(globalAlignmentSettings.get_mode()) << std::endl;
	std::cout << "Output Cycles:            ";
	for ( auto cycle : globalAlignmentSettings.get_output_cycles() ) {
		std::cout << cycle << " ";
	}
	std::cout << std::endl;
	std::cout << std::endl;
}
