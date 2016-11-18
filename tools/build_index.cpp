#include <boost/program_options.hpp>

#include "../lib/headers.h"
#include "../lib/definitions.h"
#include "../lib/kindex.h"

namespace po = boost::program_options;

std::string license =
"Copyright (c) 2015-2016, Martin S. Lindner and the HiLive contributors. See CONTRIBUTORS for more info.\n"
"All rights reserved.\n"
"\n"
"Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:\n"
"\n"
"1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.\n"
"\n"
"2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.\n"
"\n"
"3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.\n"
"\n"
"THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.";



AlignmentSettings globalAlignmentSettings; // for hard coded gapped kmer structure

int main(int argc, char* argv[]) {

  unsigned trim;
  bool do_not_convert_spaces;
  bool trim_ids;

  // setting up the command line interface
  po::options_description general("General");
  general.add_options()
    ("help,h", "Print this help message and exit")
    ("license", "Print licensing information and exit");

  po::options_description parameters("Parameters");
  parameters.add_options()
    ("INPUT", po::value<std::string>()->required(), "Input reference genome (fasta file)");

  po::options_description options("Options");
  options.add_options()
    ("outfile,o", po::value<std::string>(), "Set output file name [Default: INPUT.kix]")
    ("trim,t", po::value<unsigned>(&trim)->default_value(0), "Ignore k-mers with more than t occurrences. [Default: no limit]")
    ("do-not-convert-spaces", po::bool_switch(&do_not_convert_spaces)->default_value(false), "Do not convert all spaces in reference ids to underscores [Default: converting is on]")
    ("trim-after-space", po::bool_switch(&trim_ids)->default_value(false), "Trim all reference ids after first space [Default: false]");

  po::options_description cmdline_options;
  cmdline_options.add(general).add(parameters).add(options);

  po::options_description visible_options;
  visible_options.add(general).add(options);

  std::stringstream help_message;
  help_message << "HiLive index builder v"<< HiLive_VERSION_MAJOR << "." << HiLive_VERSION_MINOR << std::endl;
  help_message << "Index creation tool for HiLive - Realtime Alignment of Illumina Reads" << std::endl;
  help_message << "Copyright (c) 2015, Martin S. Lindner" << std::endl;
  help_message << "HiLive is open-source software. Check with --license for details." << std::endl << std::endl;
  help_message << "Fixed k-mer size: " << K_HiLive << std::endl << std::endl;
  help_message << "Usage: " << std::string(argv[0]) << " [options] INPUT" << std::endl;
  help_message << "  INPUT       Reference genomes in (multi-) FASTA format" << std::endl;

  help_message << visible_options << std::endl;

  po::positional_options_description p;
  p.add("INPUT", 1);

  po::variables_map vm;
  try {
    // parse arguments
    po::store(po::command_line_parser(argc, argv).
	      options(cmdline_options).positional(p).run(), vm);
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


  std::string fasta_name = vm["INPUT"].as<std::string>();

  std::string index_name;
  if (vm.count("outfile")) {
    // use the output file name if provided as argument
    index_name = vm["outfile"].as<std::string>();
  } else {
    // construct it from the input file name otherwise
    index_name = fasta_name + std::string(".kix");    
  }

  std::cout << "Creating index with K_HiLive=" << K_HiLive << " from file " << fasta_name << std::endl; 
  KixBuild* index = new KixBuild();
  index->add_fasta(fasta_name, !do_not_convert_spaces, trim_ids);

  if (trim > 0) {
    uint64_t trimmed = index->trim(trim);
    std::cout << "Removed " << trimmed << " k-mer positions from the database." << std::endl;
  }

  
  std::cout << "Writing index to file " << index_name << std::endl;
  index->serialize_file(index_name);

  delete index;
} 
