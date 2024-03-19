cmake_minimum_required(VERSION 3.21)
project(isoseqAligner
	VERSION 0.1
	LANGUAGES CXX
)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_EXTENSIONS OFF)
if (NOT DEFINED CMAKE_CXX_STANDARD)
	set(CMAKE_CXX_STANDARD 14)
endif()

# Adding a Profile build type
get_property(isMultiConfig GLOBAL
	PROPERTY GENERATOR_IS_MULTICONFIG
)
if(isMultiConfig)
	if(NOT "Profile" IN_LIST CMAKE_CONFIGURATION_TYPES)
		list(APPEND CMAKE_CONFIGURATION_TYPES Profile)
	endif()
	if(NOT "Test" IN_LIST CMAKE_CONFIGURATION_TYPES)
		list(APPEND CMAKE_CONFIGURATION_TYPES Test)
	endif()
else()
	set(allowedBuildTypes Debug Release Profile Test)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
		STRINGS "${allowedBuildTypes}"
	)
	if(NOT CMAKE_BUILD_TYPE)
		set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE) # default build type is Release
	elseif(NOT CMAKE_BUILD_TYPE IN_LIST allowedBuildTypes)
		message(FATAL_ERROR
			"Unknown build type ${CMAKE_BUILD_TYPE}"
		)
	endif()
endif()
#
# Set warnings
#
set(CLANG_WARNINGS
	-Wall
	-Wextra
	-Wconversion
	-Wpedantic
	-Wshadow
	-Wnon-virtual-dtor
	-Wold-style-cast
	-Wcast-align
	-Woverloaded-virtual
	-Wsign-conversion
	-Wnull-dereference
	-Wdouble-promotion
	-Wformat=2
	-Wimplicit-fallthrough
)
set(GCC_WARNINGS
	${CLANG_WARNINGS} 
	-Wmisleading-indentation
	-Wduplicated-cond
	-Wduplicated-branches
	-Wlogical-op
	-Wuseless-cast
)

if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
	set(PROJECT_WARNINGS_CXX ${CLANG_WARNINGS})
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(PROJECT_WARNINGS_CXX ${GCC_WARNINGS})
else()
	message(AUTHOR_WARNING "No compiler warnings set for CXX compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()
include(ExternalProject)
ExternalProject_Add(htslib
	#SOURCE_DIR             htslib
	GIT_REPOSITORY         https://github.com/samtools/htslib
	GIT_TAG                6a7d33abc6cae840023868ccdd946d0d8759f259
	GIT_SUBMODULES_RECURSE TRUE
	GIT_PROGRESS           TRUE
	CONFIGURE_COMMAND      autoreconf -i <SOURCE_DIR>
	COMMAND                <SOURCE_DIR>/configure
	BUILD_COMMAND          make
	INSTALL_COMMAND        ""  # no need to install system-wide
	EXCLUDE_FROM_ALL       TRUE
	SYSTEM                 # to disable warnings
)
ExternalProject_Get_Property(htslib BINARY_DIR)
set(htslib_BINARY_DIR ${BINARY_DIR})
ExternalProject_Get_Property(htslib SOURCE_DIR)
set(htslib_SOURCE_DIR ${SOURCE_DIR})

#
# Tests with Catch2
#
option(BUILD_TESTS
	"Build tests"
	OFF
)
if(PROJECT_IS_TOP_LEVEL AND BUILD_TESTS)
	Include(FetchContent)
	set (FETCHCONTENT_QUIET FALSE)
	FetchContent_Declare(
	  Catch2
	  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	  GIT_TAG        v3.5.2 # or a later release
	  GIT_PROGRESS   TRUE
	  SYSTEM                # to disable Catch2 warnings
	)
	FetchContent_MakeAvailable(Catch2)
endif()
#
# Find active available sanitizers
#
set(SANITIZER_LIST "")
if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(SANITIZERS "")
	check_cxx_compiler_flag(-fsanitize=address    HAS_ADDRESS_SANITIZER)
	check_cxx_compiler_flag(-fsanitize=leak       HAS_LEAK_SANITIZER)
	#check_cxx_compiler_flag(-fsanitize=undefined  HAS_UNDEFINED_SANITIZER)
	if(${HAS_ADDRESS_SANITIZER})
		list(APPEND SANITIZERS "address")
	endif()
	if(${HAS_LEAK_SANITIZER})
		list(APPEND SANITIZERS "leak")
	endif()
	# weird linking errors wiht UBSAN
	#if(${HAS_UNDEFINED_SANITIZER})
	#	list(APPEND SANITIZERS "undefined")
	#endif()
	list(JOIN
		SANITIZERS
		","
		SANITIZER_LIST)
endif()

set(CMAKE_CXX_FLAGS_RELEASE
	"-O3 -DNDEBUG"
)

set(CMAKE_CXX_FLAGS_PROFILE
	"-g -p -O3 -DNDEBUG"
)

set(CMAKE_CXX_FLAGS_DEBUG
	"-g -O1"
)

set(CMAKE_CXX_FLAGS_TEST
	"-g -O3"
)

include(GNUInstallDirs)
# executable
#add_executable(ldblocks
#	${CMAKE_CURRENT_LIST_DIR}/apps/ldblocks.cpp
#)
#target_include_directories(ldblocks
#	PRIVATE ${brutilities_SOURCE_DIR}/include
#	PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include
#)
#target_link_libraries(ldblocks
#	PRIVATE vash
#)
#target_compile_options(ldblocks
#	PRIVATE ${PROJECT_WARNINGS_CXX}
#)
#if(BUILD_TESTS)
#	target_compile_options(ldblocks
#		PRIVATE -fsanitize=${SANITIZER_LIST}
#	)
#endif()
#set_target_properties(ldblocks PROPERTIES
#	CXX_STANDARD_REQUIRED ON
#)
#install(TARGETS ldblocks
#	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#)

# library
#add_library(vash
#	${CMAKE_CURRENT_LIST_DIR}/src/gvarHash.cpp
#	${CMAKE_CURRENT_LIST_DIR}/src/vashFunctions.cpp
#	${CMAKE_CURRENT_LIST_DIR}/src/similarityMatrix.cpp
#	${brutilities_SOURCE_DIR}/src/random.cpp
#)
#target_include_directories(vash
#	PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include
#	PRIVATE ${brutilities_SOURCE_DIR}/include
#)
#set_target_properties(vash PROPERTIES
#	POSITION_INDEPENDENT_CODE ON
#	CXX_STANDARD_REQUIRED ON
#)
#target_compile_options(vash
#	PRIVATE ${PROJECT_WARNINGS_CXX} -mbmi -mbmi2 -mpopcnt -pthread
#)
#if(PROJECT_IS_TOP_LEVEL AND BUILD_TESTS)
#	target_compile_options(vash
#		PRIVATE -fsanitize=${SANITIZER_LIST}
#	)
#endif()
# no install for now

# tests
if(PROJECT_IS_TOP_LEVEL AND BUILD_TESTS)
	add_executable(tests
		tests/tests.cpp
	)
	add_dependencies(tests
		htslib
	)
	target_link_libraries(tests
		PRIVATE ${htslib_BINARY_DIR}/libhts.so
		PRIVATE Catch2::Catch2WithMain
	)
	target_include_directories(tests
		PRIVATE ${htslib_SOURCE_DIR}/htslib
		PRIVATE include
	)
	include(CTest)
	include(Catch)
	catch_discover_tests(tests)
endif()
