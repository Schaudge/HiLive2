#include <boost/program_options.hpp>
#include "headers.h"
#include "definitions.h"
#include "global_variables.h"
#include "parallel.h"

namespace po = boost::program_options;


/**
 * Interface for the argument parsers of the different executables.
 * @author Tobias Loka
 */
class ArgumentParser {

protected:

	/** Number of command line arguments. */
	int argc;

	/** List of command line arguments. */
	char const ** argv;

	/** Program license. */
	std::string license = "Copyright (c) 2015-2016, Martin S. Lindner and the HiLive contributors. See CONTRIBUTORS for more info.\n"
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

	/** Help output (differs for each executable). */
	std::string help;

	/**
	 * (virtual) function that is called from the executable's main() to parse the command line arguments.
	 * @return Error code ( 0: continue program; 1: exit with success; -1: exit with error)
	 * @author Tobias Loka
	 */
	virtual int parseCommandLineArguments() = 0;

	/**
	 * Print the license.
	 * @author Tobias Loka
	 */
	void printLicense(){ std::cout << license << std::endl; };

	/**
	 * Init the help text by using the options that should be visible to the user.
	 * @param visible_options The options that shall be visible for the user.
	 * @return help message as string
	 * @author Tobias Loka
	 */
	virtual void init_help(po::options_description visible_options) = 0;

	/**
	 * Print the license.
	 * @author Tobias Loka
	 */
	void printHelp(){ std::cout << help << std::endl; };

	/**
	 * Report the most important settings to the console.
	 * @author Martin Lindner
	 */
	virtual void report() = 0;

public:

	/**
	 * Default constructor.
	 * @argC Number of command line arguments.
	 * @argV List of command line arguments.
	 * @author Tobias Loka
	 */
	explicit ArgumentParser(int argC, char const ** argV);

	/**
	 * Virtual destructor.
	 * @author Tobias Loka
	 */
	virtual ~ArgumentParser(){};

};

/**
 * Class to parse arguments for HiLive build.
 */
class BuildIndexArgumentParser : public ArgumentParser {

	uint16_t kmer_weight;

	std::vector<unsigned> gap_positions = {};

	/**
	 * Use the constructor of the inherited ArgumentParser class.
	 */
	using ArgumentParser::ArgumentParser;

	/**
	 * General options of HiLive build.
	 * @return Option descriptor containing all general options that can be set by the user.
	 * @author Martin Lindner
	 */
	po::options_description general_options();

	/**
	 * Positional options of HiLive build.
	 * @return Option descriptor containing all positional options that must be set by the user.
	 * @author Martin Lindner
	 */
	po::options_description positional_options();

	/**
	 * Build options of HiLive build.
	 * @return Option descriptor containing all positional options that must be set by the user.
	 * @author Martin Lindner
	 */
	po::options_description build_options();

	/**
	 * Set all variables for the positional command arguments.
	 * @param vm The variables map containing the user parameters.
	 * @return true on success, false otherwise
	 * @author Tobias Loka
	 */
	bool set_positional_variables(po::variables_map vm);

	/**
	 * Set all variables for the build arguments.
	 * @param vm The variables map containing the user parameters.
	 * @return true on success, false otherwise
	 * @author Tobias Loka
	 */
	bool set_build_variables(po::variables_map vm);

	void report() override;

	void init_help(po::options_description visible_options) override;

public:

	// name of the index file
	std::string index_name;

	// name of the input fasta file
	std::string fasta_name;

	// trimming parameter
	unsigned trim;

	// do_not_convert_spaces_switch
	bool do_not_convert_spaces;

	// trim_ids switch
	bool trim_ids;

	int parseCommandLineArguments() override;

};


/**
 * Class to parse arguments for HiLive.
 */
class HiLiveArgumentParser : public ArgumentParser {

	/**
	 * Use the constructor of the inherited ArgumentParser class.
	 */
	using ArgumentParser::ArgumentParser;

	/**
	 * General options of HiLive.
	 * @return Option descriptor containing all general options that can be set by the user.
	 * @author Martin Lindner
	 */
	po::options_description general_options();

	/**
	 * Positional options of HiLive.
	 * @return Option descriptor containing all positional options that must be set by the user.
	 * @author Martin Lindner
	 */
	po::options_description positional_options();

	/**
	 * I/O options of HiLive.
	 * @return Option descriptor containing all I/O options that can be set by the user.
	 * @author Martin Lindner
	 */
	po::options_description io_options();

	/**
	 * Alignment options of HiLive.
	 * @return Option descriptor containing all alignment options that can be set by the user.
	 * @author Martin Lindner
	 */
	po::options_description alignment_options();

	/**
	 * Technical options of HiLive.
	 * @return Option descriptor containing all technical options that can be set by the user.
	 * @author Martin Lindner
	 */
	po::options_description technical_options();

	/**
	 * Set all variables for the positional command arguments.
	 * @param vm The variables map containing the user parameters.
	 * @return true on success, false otherwise
	 * @author Tobias Loka
	 */
	bool set_positional_variables(po::variables_map vm);

	/**
	 * Set all variables for the I/O settings.
	 * @param vm The variables map containing the user parameters.
	 * @return true on success, false otherwise
	 * @author Tobias Loka
	 */
	bool set_io_variables(po::variables_map vm);

	/**
	 * Set all variables for the alignment settings.
	 * @param vm The variables map containing the user parameters.
	 * @return true on success, false otherwise
	 * @author Tobias Loka
	 */
	bool set_alignment_variables(po::variables_map vm);

	/**
	 * Set all variables for the technical settings.
	 * @param vm The variables map containing the user parameters.
	 * @return true on success, false otherwise
	 * @author Tobias Loka
	 */
	bool set_technical_variables(po::variables_map vm);

	/**
	 * Check all paths that are relevant for the functionality of HiLive.
	 * @return true if all paths and files are accessible
	 * @author Jakob Schulze
	 */
	bool checkPaths();

	/**
	 * Parse the user parameter to define the read fragments (--reads,-r).
	 * @param readsArg Vector of strings containing the read fragments
	 * @return true on success, false otherwise
	 * @author Tobias Loka
	 */
	bool parseReadsArgument(std::vector < std::string > readsArg);

	/**
	 * Parse the user parameter to define the barcodes (--barcodes,-b).
	 * @param barcodeArg Vector of strings containing the barcodes. Barcode components are delimited by a "-".
	 * @return true on success, false otherwise
	 * @author Tobias Loka
	 */
	bool parseBarcodeArgument(std::vector < std::string > barcodeArg);

	/**
	 * Parse Lanes, Tiles and read fragments from a RunInfo.xml file.
	 * @param vm The variables map containing the user parameters.
	 * @return true on success, false otherwise
	 * @author Jakob Schulze
	 */
	bool parseRunInfo(po::variables_map vm);

	void report() override;

	void init_help(po::options_description visible_options) override;

public:

	int parseCommandLineArguments() override;

};
