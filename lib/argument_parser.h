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

	/** Vector of required options. */
	std::vector<std::string> required_options;

	/** Map of command line arguments. */
	po::variables_map cmd_settings;

	/** Settings property tree obtained from the runInfo. */
	boost::property_tree::ptree runInfo_settings;

	/** Settings property tree obtained from an settings input file. */
	boost::property_tree::ptree input_settings;

	/** Bool describing whether an input settings file was specified. */
	bool has_input_settings = false;

	/** Program license. */
	std::string license = "Copyright (c) 2015-2017, Martin S. Lindner and the HiLive contributors. See CONTRIBUTORS for more info.\n"
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

	virtual void set_required_parameters() {};

	/**
	 * Check if an option is required.
	 * @param vm_key Key for the virtual map for the parameter to check (equals cmd line option name).
	 * @return true, if the option is required.
	 */
	bool isRequired( std::string vm_key) {
		if ( std::find(required_options.begin(), required_options.end(), vm_key) == required_options.end() )
			return false;
		return true;
	}

	/**
	 * Set an option in the globalAlignmentSettings.
	 * Thereby, the different input sources have the following priority:
	 * 	1. Command line argument
	 * 	2. RunInfo file
	 * 	3. Input settings file
	 * @param vm_key Key for the command line's variables map.
	 * @param settings_key Key for the RunInfo and Input settings property tree (must be similar for both files)
	 * @param default_value A default value that is set if the variable is set in none of the input sources (this is not considered if "required" is true)
	 * @param function Function that is called to set the variable in globalAlignmentSettings. Must be of type void &AlignmentSettings::*(T).
	 * @param required true, if the option must be set by the user (from one of the input sources)
	 * @author Tobias Loka
	 */
	template<class T> void set_option(std::string vm_key, std::string settings_key, T default_value, void (AlignmentSettings::*function)(T)) {
		set_option_impl(vm_key, settings_key, default_value, function, isRequired(vm_key), static_cast<T*>(0));
	}

	/**
	 * General implementation of set_option(...).
	 */
	template<class T> void set_option_impl(std::string vm_key, std::string settings_key, T default_value, void (AlignmentSettings::*function)(T), bool required, T*) {

		T value = default_value;
		bool was_set = false;

		boost::optional<T> isv = input_settings.get_optional<T>(settings_key);
		boost::optional<T> rsv = runInfo_settings.get_optional<T>(settings_key);


		// User parameter -> first priority
		if ( cmd_settings.count(vm_key) ) {
			value = cmd_settings[vm_key].as<T>();
			was_set = true;
		}

		// Settings file -> second priority
		else if ( isv ) {
			value = isv.get();
			was_set = true;
		}

		// RunInfo file -> third priority
		else if ( rsv ) {
			value = rsv.get();
			was_set = true;
		}

		// Throw exception if unset
		if ( required && !was_set )
			throw po::required_option(vm_key);

		// Otherwise set value
		auto binded_function = std::bind(function, &globalAlignmentSettings, std::placeholders::_1);
		binded_function(value);

	}

	/**
	 * Overload of set_option_impl for bool data type.
	 */
	void set_option_impl(std::string vm_key, std::string settings_key, bool default_value, void (AlignmentSettings::*function)(bool), bool required, bool *) {

		bool value = default_value;
		bool was_set = false;

		boost::optional<bool> isv = input_settings.get_optional<bool>(settings_key);
		boost::optional<bool> rsv = runInfo_settings.get_optional<bool>(settings_key);

		if ( cmd_settings[vm_key].as<bool>() != default_value ) {
			value = !default_value;
			was_set = true;
		}

		else if ( isv && isv.get() != default_value ) {
			value = !default_value;
			was_set = true;
		}

		else if ( rsv && rsv.get() != default_value ) {
			value = !default_value;
			was_set = true;
		}


		if ( required && !was_set )
			throw po::required_option(vm_key);

		// Otherwise set value
		auto binded_function = std::bind(function, &globalAlignmentSettings, std::placeholders::_1);
		binded_function(value);

	}

	/** Overload of set_option_impl for std::vector data types. */
	template<class T> void set_option_impl(std::string vm_key, std::string settings_key, std::vector<T> default_value, void (AlignmentSettings::*function)(std::vector<T>), bool required, std::vector<T> *) {

		std::vector<T> value;
		bool was_set = false;

		auto sub_isv = input_settings.get_child_optional(settings_key);
		auto sub_rsv = input_settings.get_child_optional(settings_key);


		// User parameter -> first priority
		if ( cmd_settings.count(vm_key) ) {
			value = cmd_settings[vm_key].as<std::vector<T>>();
			was_set = true;
		}

		// Settings file -> third priority
		else if ( sub_isv && sub_isv.get().count("el") ) {
			for ( auto& v : sub_isv.get() ) {
				if ( v.first == "el" )
					value.push_back(v.second.get_value<T>());
			}
			was_set = true;
		}

		// RunInfo file -> second priority
		else if ( sub_rsv && sub_rsv.get().count("el") ) {
			for ( auto& v : sub_rsv.get() ) {
				if ( v.first == "el" )
					value.push_back(v.second.get_value<T>());
			}
			was_set = true;
		}

		// Throw exception if unset
		if ( required && !was_set )
			throw po::required_option(vm_key);

		if ( !was_set) {
			value = default_value;
		}

		// Otherwise set value
		auto binded_function = std::bind(function, &globalAlignmentSettings, std::placeholders::_1);
		binded_function(value);
	}

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

protected:

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
	 * Scorings scheme of HiLive.
	 * @return Option descriptor containing all scoring options that can be set by the user.
	 * @author Martin Lindner
	 */
	po::options_description scoring_options();

	/**
	 * Technical options of HiLive.
	 * @return Option descriptor containing all technical options that can be set by the user.
	 * @author Martin Lindner
	 */
	po::options_description technical_options();

	/**
	 * Check all paths that are relevant for the functionality of HiLive.
	 * @return true if all paths and files are accessible
	 * @author Jakob Schulze
	 */
	bool checkPaths();

	/**
	 * Parse Lanes, Tiles and read fragments from a RunInfo.xml file.
	 * @param vm The variables map containing the user parameters.
	 * @return true on success, false otherwise
	 * @author Jakob Schulze
	 */
	bool parseRunInfo(po::variables_map vm);

	virtual void report() override;

	void init_help(po::options_description visible_options) override;

	bool set_options();

	virtual void set_required_parameters() override { required_options = {"BC_DIR", "INDEX", "CYCLES"}; }

public:

	int parseCommandLineArguments() override;

};

/**
 * Class to parse arguments for HiLive out.
 */

class HiLiveOutArgumentParser : public HiLiveArgumentParser {

	using HiLiveArgumentParser::HiLiveArgumentParser;

	void init_help(po::options_description visible_options) override;

	void report() override;

	void set_required_parameters() override { required_options = {"settings", "INDEX"}; };
};

//class HiLiveOutArgumentParser : public ArgumentParser {
//
//	/**
//	 * Use the constructor of the inherited ArgumentParser class.
//	 */
//	using ArgumentParser::ArgumentParser;
//
//	/**
//	 * General options of HiLive build.
//	 * @return Option descriptor containing all general options that can be set by the user.
//	 * @author Martin Lindner
//	 */
//	po::options_description general_options();
//
//	/**
//	 * Positional options of HiLive build.
//	 * @return Option descriptor containing all positional options that must be set by the user.
//	 * @author Martin Lindner
//	 */
//	po::options_description positional_options();
//
//	/**
//	 * Build options of HiLive build.
//	 * @return Option descriptor containing all positional options that must be set by the user.
//	 * @author Martin Lindner
//	 */
//	po::options_description output_options();
//
//	/**
//	 * Set all variables for the positional command arguments.
//	 * @param vm The variables map containing the user parameters.
//	 * @return true on success, false otherwise
//	 * @author Tobias Loka
//	 */
//	bool set_positional_variables(po::variables_map vm);
//
//	/**
//	 * Set all variables for the build arguments.
//	 * @param vm The variables map containing the user parameters.
//	 * @return true on success, false otherwise
//	 * @author Tobias Loka
//	 */
//	bool set_output_variables(po::variables_map vm);
//
//	void report() override;
//
//	void init_help(po::options_description visible_options) override;
//
//	bool set_options() override;
//
//public:
//
//	int parseCommandLineArguments() override;
//
//};
