cmake_minimum_required(VERSION 3.21)
project(isoSeqQC
	VERSION 0.2
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
# executables
add_executable(exoncoverage
	${CMAKE_CURRENT_LIST_DIR}/apps/exoncoverage.cpp
)
target_include_directories(exoncoverage
	PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include
)
target_include_directories(exoncoverage SYSTEM
	PRIVATE ${htslib_SOURCE_DIR}/htslib
)
target_link_libraries(exoncoverage
	PRIVATE isoseqAln
)
target_compile_options(exoncoverage
	PRIVATE ${PROJECT_WARNINGS_CXX}
)
if(BUILD_TESTS)
	target_compile_options(exoncoverage
		PRIVATE -fsanitize=${SANITIZER_LIST}
	)
endif()
set_target_properties(exoncoverage PROPERTIES
	CXX_STANDARD_REQUIRED ON
)
install(TARGETS exoncoverage
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

add_executable(partialmaps
	${CMAKE_CURRENT_LIST_DIR}/apps/partialmaps.cpp
)
target_include_directories(partialmaps
	PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include
)
target_include_directories(partialmaps SYSTEM
	PRIVATE ${htslib_SOURCE_DIR}/htslib
)
target_link_libraries(partialmaps
	PRIVATE isoseqAln
)
target_compile_options(partialmaps
	PRIVATE ${PROJECT_WARNINGS_CXX}
)
if(BUILD_TESTS)
	target_compile_options(partialmaps
		PRIVATE -fsanitize=${SANITIZER_LIST}
	)
endif()
set_target_properties(partialmaps PROPERTIES
	CXX_STANDARD_REQUIRED ON
)
install(TARGETS partialmaps
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

add_executable(addremaps
	${CMAKE_CURRENT_LIST_DIR}/apps/addremaps.cpp
)
target_include_directories(addremaps
	PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include
)
target_include_directories(addremaps SYSTEM
	PRIVATE ${htslib_SOURCE_DIR}/htslib
)
target_link_libraries(addremaps
	PRIVATE isoseqAln
)
target_compile_options(addremaps
	PRIVATE ${PROJECT_WARNINGS_CXX}
)
if(BUILD_TESTS)
	target_compile_options(addremaps
		PRIVATE -fsanitize=${SANITIZER_LIST}
	)
endif()
set_target_properties(addremaps PROPERTIES
	CXX_STANDARD_REQUIRED ON
)
install(TARGETS addremaps
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
# library
add_library(isoseqAln
	${CMAKE_CURRENT_LIST_DIR}/src/isoseqAlgn.cpp
	${CMAKE_CURRENT_LIST_DIR}/src/helperFunctions.cpp
)
add_dependencies(isoseqAln htslib)
target_link_libraries(isoseqAln
	PRIVATE ${htslib_BINARY_DIR}/libhts.so
)
target_include_directories(isoseqAln
	PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include
)
target_include_directories(isoseqAln SYSTEM
	PRIVATE ${htslib_SOURCE_DIR}/htslib
)
set_target_properties(isoseqAln PROPERTIES
	POSITION_INDEPENDENT_CODE ON
	CXX_STANDARD_REQUIRED ON
)
if(PROJECT_IS_TOP_LEVEL AND BUILD_TESTS)
	target_compile_options(isoseqAln
		PRIVATE -fsanitize=${SANITIZER_LIST}
	)
endif()
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
		PRIVATE isoseqAln
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
