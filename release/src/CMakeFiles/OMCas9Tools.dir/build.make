# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/cshi/tools/OMwithCas9_Tools

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cshi/tools/OMwithCas9_Tools/release

# Include any dependencies generated for this target.
include src/CMakeFiles/OMCas9Tools.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/OMCas9Tools.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/OMCas9Tools.dir/flags.make

src/CMakeFiles/OMCas9Tools.dir/main.cpp.o: src/CMakeFiles/OMCas9Tools.dir/flags.make
src/CMakeFiles/OMCas9Tools.dir/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cshi/tools/OMwithCas9_Tools/release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/OMCas9Tools.dir/main.cpp.o"
	cd /home/cshi/tools/OMwithCas9_Tools/release/src && /home/cshi/tools/Miniconda3/envs/OM/bin/x86_64-conda-linux-gnu-c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/OMCas9Tools.dir/main.cpp.o -c /home/cshi/tools/OMwithCas9_Tools/src/main.cpp

src/CMakeFiles/OMCas9Tools.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/OMCas9Tools.dir/main.cpp.i"
	cd /home/cshi/tools/OMwithCas9_Tools/release/src && /home/cshi/tools/Miniconda3/envs/OM/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cshi/tools/OMwithCas9_Tools/src/main.cpp > CMakeFiles/OMCas9Tools.dir/main.cpp.i

src/CMakeFiles/OMCas9Tools.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/OMCas9Tools.dir/main.cpp.s"
	cd /home/cshi/tools/OMwithCas9_Tools/release/src && /home/cshi/tools/Miniconda3/envs/OM/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cshi/tools/OMwithCas9_Tools/src/main.cpp -o CMakeFiles/OMCas9Tools.dir/main.cpp.s

# Object files for target OMCas9Tools
OMCas9Tools_OBJECTS = \
"CMakeFiles/OMCas9Tools.dir/main.cpp.o"

# External object files for target OMCas9Tools
OMCas9Tools_EXTERNAL_OBJECTS =

src/OMCas9Tools: src/CMakeFiles/OMCas9Tools.dir/main.cpp.o
src/OMCas9Tools: src/CMakeFiles/OMCas9Tools.dir/build.make
src/OMCas9Tools: src/liblib_core.a
src/OMCas9Tools: src/CMakeFiles/OMCas9Tools.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/cshi/tools/OMwithCas9_Tools/release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable OMCas9Tools"
	cd /home/cshi/tools/OMwithCas9_Tools/release/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/OMCas9Tools.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/OMCas9Tools.dir/build: src/OMCas9Tools

.PHONY : src/CMakeFiles/OMCas9Tools.dir/build

src/CMakeFiles/OMCas9Tools.dir/clean:
	cd /home/cshi/tools/OMwithCas9_Tools/release/src && $(CMAKE_COMMAND) -P CMakeFiles/OMCas9Tools.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/OMCas9Tools.dir/clean

src/CMakeFiles/OMCas9Tools.dir/depend:
	cd /home/cshi/tools/OMwithCas9_Tools/release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cshi/tools/OMwithCas9_Tools /home/cshi/tools/OMwithCas9_Tools/src /home/cshi/tools/OMwithCas9_Tools/release /home/cshi/tools/OMwithCas9_Tools/release/src /home/cshi/tools/OMwithCas9_Tools/release/src/CMakeFiles/OMCas9Tools.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/OMCas9Tools.dir/depend

