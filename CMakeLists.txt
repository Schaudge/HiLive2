# Header
cmake_minimum_required (VERSION 2.8)
project (HiLive)
include(ExternalProject)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -pthread -W -Wall -std=gnu++11 -O3")

# Set the k-mer length
set (HiLive_K 15 CACHE INT "Set the k-mer length for index and mapper")

# Set the version number 
set (HiLive_VERSION_MAJOR 0)
set (HiLive_VERSION_MINOR 2)

configure_file (
  "${PROJECT_SOURCE_DIR}/lib/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"  )

# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")



include_directories(lib)

set(LIBS tools
         alnread
	 alnstream
	 illumina_parsers
	 kindex
	 parallel
	 bamout)

set(LIB_FILES "")
foreach (lib ${LIBS})
    list(APPEND LIB_FILES "lib/${lib}.cpp")
endforeach()

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME ON) 
FIND_PACKAGE( Boost COMPONENTS system filesystem program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

find_package( ZLIB REQUIRED )

# manually set the path for lz4 library
set (LZ4_PATH /usr/local/lib CACHE PATH "Manually set the path to the lz4 library")
LINK_DIRECTORIES(${LZ4_PATH})

# download & build the Samtools library
file(DOWNLOAD "http://downloads.sourceforge.net/project/samtools/samtools/1.2/samtools-1.2.tar.bz2" ${CMAKE_CURRENT_BINARY_DIR}/samtools-1.2.tar.bz2)
execute_process(COMMAND tar xf ${CMAKE_CURRENT_BINARY_DIR}/samtools-1.2.tar.bz2 )
execute_process(COMMAND make -s -C samtools-1.2 )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR}/samtools-1.2 )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR}/samtools-1.2/htslib-1.2.1 )
LINK_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR}/samtools-1.2 ${CMAKE_CURRENT_BINARY_DIR}/samtools-1.2/htslib-1.2.1 )

# create the executables
add_library(hilivelib ${LIB_FILES})

add_executable(hilive tools/hilive.cpp )
target_link_libraries(hilive hilivelib ${ZLIB_LIBRARIES} ${Boost_LIBRARIES} lz4 hts)

add_executable(hilive-build tools/build_index.cpp )
target_link_libraries(hilive-build hilivelib ${ZLIB_LIBRARIES} ${Boost_LIBRARIES} lz4)
