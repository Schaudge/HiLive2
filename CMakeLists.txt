##############
### Header ###

cmake_minimum_required (VERSION 2.8)
project (HiLive)

# Set the version number 
add_definitions(-DHiLive_VERSION_MAJOR=1)
add_definitions(-DHiLive_VERSION_MINOR=0)

# Set flags for compilation
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -g -pthread -W -Wall -std=gnu++14 -O0")

# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")


#############################
### setup Boost libraries ###

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME ON) 
find_package( Boost COMPONENTS system filesystem program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )


############################
### setup Zlib and lz4 libraries ###

find_package( ZLIB REQUIRED )
#set (LZ4_PATH /usr/local/lib) 		# possibly adjust this to [pathToLz4]/lib if using downloaded lz4 source code
include_directories(${LZ4_PATH})
link_directories(${LZ4_PATH})
set(CompressionLibs "${ZLIB_LIBRARIES};lz4")


#############################
### setup seqan libraries ### needs to be done AFTER searching for Zlib

set (CMAKE_MODULE_PATH "/usr/local/lib/seqan/util/cmake") # adjust this to [pathToSeqanCode]/util/cmake
set (SEQAN_INCLUDE_PATH "/usr/local/lib/seqan/include/") # adjust this to [pathToSeqanCode]/include

# Configure SeqAn, enabling features for libbz2 and zlib.
#set (SEQAN_FIND_DEPENDENCIES ZLIB BZip2) # original version from seqan tutorial
set (SEQAN_FIND_DEPENDENCIES BZip2)
find_package (SeqAn REQUIRED)

# Add include directories, defines, and flags for SeqAn (and its dependencies).
include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")


##############################
### setup HiLive libraries ###
include_directories("${PROJECT_SOURCE_DIR}/lib")
# make a list of HiLive libraries
set (LIB_NAMES tools_static tools alnread alnstream illumina_parsers kindex parallel argument_parser)
set(LIB_LIST "")
foreach (x ${LIB_NAMES})
	list(APPEND LIB_LIST "lib/${x}.cpp")
endforeach()
add_library(HiLiveLibs ${LIB_LIST})


#############################
### Build the executables ###

add_executable (hilive tools/hilive.cpp)
target_link_libraries (hilive HiLiveLibs ${CompressionLibs} ${Boost_LIBRARIES} ${SEQAN_LIBRARIES})

add_executable(hilive-build tools/build_index.cpp )
target_link_libraries(hilive-build HiLiveLibs ${CompressionLibs} ${Boost_LIBRARIES} ${SEQAN_LIBRARIES})

add_executable(hilive-out tools/hilive_out.cpp )
target_link_libraries(hilive-out HiLiveLibs ${CompressionLibs} ${Boost_LIBRARIES} ${SEQAN_LIBRARIES})


#####################
### for debugging ###

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
	#message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
