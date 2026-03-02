find_package(Doxygen REQUIRED)
find_package(LATEX
	COMPONENTS PDFLATEX
	REQUIRED
)

message("pdflatex: ${PDFLATEX_COMPILER}")

FetchContent_Declare(
	awesome
	GIT_REPOSITORY https://github.com/jothepro/doxygen-awesome-css.git
	GIT_TAG        v2.3.4 # or a later release
	GIT_PROGRESS   TRUE
)
FetchContent_MakeAvailable(awesome)

file(GLOB PUBLIC_HEADERS ${PROJECT_SOURCE_DIR}/include/*.hpp)

set(DOXYGEN_INPUT_DIR  ${PROJECT_SOURCE_DIR}/include)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/doxygen)
set(README_FILE        ${PROJECT_SOURCE_DIR}/README.md)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
set(LATEX_REFMAN_TEX   ${DOXYGEN_OUTPUT_DIR}/latex/refman.tex)
set(LATEX_REFMAN_PDF   ${DOXYGEN_OUTPUT_DIR}/latex/refman.pdf)
set(DOXYFILE_IN        ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT       ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
set(AWESOME_CSS        ${awesome_SOURCE_DIR}/doxygen-awesome.css)

# Replace variables in Doxyfile
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})
add_custom_command(
	OUTPUT            ${DOXYGEN_INDEX_FILE}
	DEPENDS           ${PUBLIC_HEADERS}
	                  ${README_FILE}
	COMMAND           ${DOXYGEN_EXECUTABLE}
	COMMAND           make --directory=${DOXYGEN_OUTPUT_DIR}/latex pdf
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	MAIN_DEPENDENCY   ${DOXYFILE_IN}
	COMMENT           "Generating docs"
	VERBATIM
)
add_custom_target(Doxygen
	ALL
	DEPENDS ${DOXYGEN_INDEX_FILE}
)
