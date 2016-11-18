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
        //("bam,B", po::bool_switch()->default_value(false), "Create BAM files instead of SAM files [Default: false]")
        ("keep-files,k", po::bool_switch()->default_value(false), "Keep intermediate alignment files [Default: false]")
        ("lanes,l", po::value< std::vector<uint16_t> >()->multitoken()->composing(), "Select lane [Default: all lanes]")
        ("tiles,t", po::value< std::vector<uint16_t> >()->multitoken()->composing(), "Select tile numbers [Default: all tiles]")
        ("barcodes,b", po::value< std::vector<std::string> >()->multitoken()->composing(), "Enumerate barcodes (must have same length) for demultiplexing, i.e. -b AGGATC -b CCCTTT [Default: no demultiplexing]");

    po::options_description alignment("Alignment settings");
    alignment.add_options()
        ("min-errors,e", po::value<CountType>()->default_value(2), "Number of errors tolerated in read alignment [Default: 2]")
        ("all-best-hit,H", po::bool_switch()->default_value(false), "Report all of the best alignments for each read")
        ("any-best-hit", po::bool_switch(), "Report one of the best alignments for each read (default)")
        ("all-best-n-scores,N", po::value<CountType>(), "Report all alignments of the N best alignment scores for each read")
        ("all-hits,A", po::bool_switch()->default_value(false), "Report all valid alignments for each read")
        ("disable-ohw-filter", po::bool_switch()->default_value(true), "disable the One-Hit Wonder filter [Default: false]")
        ("start-ohw", po::value<CountType>()->default_value(K_HiLive+5), "First cycle to apply One-Hit Wonder filter [Default: K+5]")
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


    globalAlignmentSettings.root = vm["BC_DIR"].as<std::string>();
    globalAlignmentSettings.index_fname = vm["INDEX"].as<std::string>();
    globalAlignmentSettings.rlen = vm["CYCLES"].as<CountType>();
    globalAlignmentSettings.temp_dir = vm["temp"].as<std::string>();
    //globalAlignmentSettings.write_bam = vm["bam"].as<bool>();
    globalAlignmentSettings.keep_aln_files = vm["keep-files"].as<bool>();
    globalAlignmentSettings.min_errors = vm["min-errors"].as<CountType>();

    globalAlignmentSettings.discard_ohw = vm["disable-ohw-filter"].as<bool>();
    globalAlignmentSettings.start_ohw = vm["start-ohw"].as<CountType>();
    globalAlignmentSettings.window = vm["window"].as<DiffType>();
    globalAlignmentSettings.min_qual = vm["min-quality"].as<CountType>();
    globalAlignmentSettings.compression_format = vm["compression"].as<uint8_t>();


    if (vm.count("OUTDIR"))
        globalAlignmentSettings.out_dir = vm["OUTDIR"].as<std::string>();
    else {
        if (globalAlignmentSettings.temp_dir == "") 
            globalAlignmentSettings.out_dir = globalAlignmentSettings.root;
        else
            globalAlignmentSettings.out_dir = globalAlignmentSettings.temp_dir;
    }

    if (vm.count("lanes"))
        globalAlignmentSettings.lanes = vm["lanes"].as< std::vector<uint16_t> >();
    else
        globalAlignmentSettings.lanes = all_lanes();

    if (vm.count("tiles"))
        globalAlignmentSettings.tiles = vm["tiles"].as< std::vector<uint16_t> >();
    else
        globalAlignmentSettings.tiles = all_tiles();

    globalAlignmentSettings.barcodeVector.clear();
    globalAlignmentSettings.seqlen = globalAlignmentSettings.rlen;
    if (vm.count("barcodes")) {
        globalAlignmentSettings.barcodeVector = vm["barcodes"].as< std::vector<std::string> >();
        globalAlignmentSettings.seqlen = globalAlignmentSettings.rlen - globalAlignmentSettings.barcodeVector[0].size();
    }

    if (vm["all-hits"].as<bool>()) {
        // all hits: disable other modes
        globalAlignmentSettings.any_best_hit_mode = false;
        globalAlignmentSettings.all_best_hit_mode = false;
        globalAlignmentSettings.all_best_n_scores_mode = false;
    } else if (vm["all-best-hit"].as<bool>()) {
        // enable all-best-hit mode and disable others
        globalAlignmentSettings.any_best_hit_mode = false;
        globalAlignmentSettings.all_best_hit_mode = true;
        globalAlignmentSettings.all_best_n_scores_mode = false;
    } else if (vm.count("all-best-n-scores")) {
        // enable any-best-n mode and get parameter
        globalAlignmentSettings.any_best_hit_mode = false;
        globalAlignmentSettings.all_best_hit_mode = false;
        globalAlignmentSettings.all_best_n_scores_mode = true;
        globalAlignmentSettings.best_n = vm["all-best-n-scores"].as<CountType>();
    } else { // the default behaviour
        // enable any-best-hit mode and disable others
        globalAlignmentSettings.any_best_hit_mode = true;
        globalAlignmentSettings.all_best_hit_mode = false;
        globalAlignmentSettings.all_best_n_scores_mode = false;
    }
        
    if (vm["-M"].as<bool>())
        globalAlignmentSettings.block_size = vm["block-size"].as<uint64_t>()*1024*1024;
    else if (vm["-K"].as<bool>())
        globalAlignmentSettings.block_size = vm["block-size"].as<uint64_t>()*1024;
    else
        globalAlignmentSettings.block_size = vm["block-size"].as<uint64_t>();

    if (vm.count("num-threads")) 
        globalAlignmentSettings.num_threads = vm["num-threads"].as<int>();
    else { // the default case, meaning as much as physically useful
        uint32_t n_cpu = std::thread::hardware_concurrency();
        if (n_cpu > 1)
            globalAlignmentSettings.num_threads = n_cpu;
        else
            globalAlignmentSettings.num_threads = 1;
    }


    // check paths and file names
    if (!file_exists(globalAlignmentSettings.index_fname)){
        std::cerr << "Input error: Could not find k-mer index file " << globalAlignmentSettings.index_fname << std::endl;
        return -1;
    }

    std::size_t found = globalAlignmentSettings.root.find("BaseCalls");
    if (!(found != std::string::npos && found >= globalAlignmentSettings.root.size()-10)) {
        std::cerr << "Warning: BaseCalls directory seems to be invalid: " << globalAlignmentSettings.root << std::endl;
    } 

    if (!is_directory(globalAlignmentSettings.root)){
        std::cerr << "Input error: Could not find BaseCalls directory " << globalAlignmentSettings.root << std::endl;
        return -1;
    }

    for ( uint16_t ln : globalAlignmentSettings.lanes ) {
        std::string ln_dir = globalAlignmentSettings.root;
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
    std::cout << "Running HiLive with       " << globalAlignmentSettings.num_threads << " thread(s)." << std::endl;
    std::cout << "BaseCalls directory:      " << globalAlignmentSettings.root << std::endl;
    if (globalAlignmentSettings.temp_dir != "") {
        std::cout << "Temporary directory:      " << globalAlignmentSettings.temp_dir << std::endl;
    }
    //if (!globalAlignmentSettings.write_bam)
    std::cout << "SAM output directory:     " << globalAlignmentSettings.out_dir << std::endl;
    //else
        //std::cout << "BAM output directory:     " << globalAlignmentSettings.out_dir << std::endl;
    std::cout << "Lanes:                    ";
    for ( uint16_t ln : globalAlignmentSettings.lanes )
        std::cout << ln << " ";
    std::cout << std::endl;
    std::cout << "K-mer index:              " << globalAlignmentSettings.index_fname << std::endl;
    std::cout << "Read length:              " << globalAlignmentSettings.rlen << std::endl;
    std::cout << "Mapping error:            " << globalAlignmentSettings.min_errors << std::endl;
    if (globalAlignmentSettings.any_best_hit_mode) 
        std::cout << "Mapping mode:             Any-Best-Hit-Mode" << std::endl;
    else if (globalAlignmentSettings.all_best_hit_mode) 
        std::cout << "Mapping mode:             All-Best-Hit-Mode" << std::endl;
    else if (globalAlignmentSettings.all_best_n_scores_mode) 
        std::cout << "Mapping mode:             All-Best-N-Scores-Mode with N=" << globalAlignmentSettings.best_n << std::endl;
    else
        std::cout << "Mapping mode:             All-Hits-Mode" << std::endl;
    std::cout << std::endl;

    return 0;
}
