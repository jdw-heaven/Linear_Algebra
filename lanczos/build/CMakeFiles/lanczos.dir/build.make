# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/heaven/Desktop/doc/Linear_Algebra/lanczos

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/heaven/Desktop/doc/Linear_Algebra/lanczos/build

# Include any dependencies generated for this target.
include CMakeFiles/lanczos.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/lanczos.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/lanczos.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/lanczos.dir/flags.make

CMakeFiles/lanczos.dir/lanczos.cpp.o: CMakeFiles/lanczos.dir/flags.make
CMakeFiles/lanczos.dir/lanczos.cpp.o: /home/heaven/Desktop/doc/Linear_Algebra/lanczos/lanczos.cpp
CMakeFiles/lanczos.dir/lanczos.cpp.o: CMakeFiles/lanczos.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/heaven/Desktop/doc/Linear_Algebra/lanczos/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/lanczos.dir/lanczos.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lanczos.dir/lanczos.cpp.o -MF CMakeFiles/lanczos.dir/lanczos.cpp.o.d -o CMakeFiles/lanczos.dir/lanczos.cpp.o -c /home/heaven/Desktop/doc/Linear_Algebra/lanczos/lanczos.cpp

CMakeFiles/lanczos.dir/lanczos.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lanczos.dir/lanczos.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/heaven/Desktop/doc/Linear_Algebra/lanczos/lanczos.cpp > CMakeFiles/lanczos.dir/lanczos.cpp.i

CMakeFiles/lanczos.dir/lanczos.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lanczos.dir/lanczos.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/heaven/Desktop/doc/Linear_Algebra/lanczos/lanczos.cpp -o CMakeFiles/lanczos.dir/lanczos.cpp.s

# Object files for target lanczos
lanczos_OBJECTS = \
"CMakeFiles/lanczos.dir/lanczos.cpp.o"

# External object files for target lanczos
lanczos_EXTERNAL_OBJECTS =

lanczos: CMakeFiles/lanczos.dir/lanczos.cpp.o
lanczos: CMakeFiles/lanczos.dir/build.make
lanczos: /home/heaven/Desktop/doc/libtorch/lib/libtorch.so
lanczos: /home/heaven/Desktop/doc/libtorch/lib/libc10.so
lanczos: /home/heaven/Desktop/doc/libtorch/lib/libkineto.a
lanczos: /home/heaven/Desktop/doc/libtorch/lib/libc10.so
lanczos: CMakeFiles/lanczos.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/heaven/Desktop/doc/Linear_Algebra/lanczos/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable lanczos"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lanczos.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/lanczos.dir/build: lanczos
.PHONY : CMakeFiles/lanczos.dir/build

CMakeFiles/lanczos.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/lanczos.dir/cmake_clean.cmake
.PHONY : CMakeFiles/lanczos.dir/clean

CMakeFiles/lanczos.dir/depend:
	cd /home/heaven/Desktop/doc/Linear_Algebra/lanczos/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/heaven/Desktop/doc/Linear_Algebra/lanczos /home/heaven/Desktop/doc/Linear_Algebra/lanczos /home/heaven/Desktop/doc/Linear_Algebra/lanczos/build /home/heaven/Desktop/doc/Linear_Algebra/lanczos/build /home/heaven/Desktop/doc/Linear_Algebra/lanczos/build/CMakeFiles/lanczos.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/lanczos.dir/depend

