# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_SOURCE_DIR = /home/mottierr/Codes/Gar6more/Gar6more2d

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mottierr/Codes/Gar6more/Gar6more2d/build

# Include any dependencies generated for this target.
include CMakeFiles/Gar6more2D.out.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/Gar6more2D.out.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Gar6more2D.out.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Gar6more2D.out.dir/flags.make

CMakeFiles/Gar6more2D.out.dir/lib/bin/Gar6more2D.F90.o: CMakeFiles/Gar6more2D.out.dir/flags.make
CMakeFiles/Gar6more2D.out.dir/lib/bin/Gar6more2D.F90.o: ../lib/bin/Gar6more2D.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mottierr/Codes/Gar6more/Gar6more2d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/Gar6more2D.out.dir/lib/bin/Gar6more2D.F90.o"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/mottierr/Codes/Gar6more/Gar6more2d/lib/bin/Gar6more2D.F90 -o CMakeFiles/Gar6more2D.out.dir/lib/bin/Gar6more2D.F90.o

CMakeFiles/Gar6more2D.out.dir/lib/bin/Gar6more2D.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/Gar6more2D.out.dir/lib/bin/Gar6more2D.F90.i"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/mottierr/Codes/Gar6more/Gar6more2d/lib/bin/Gar6more2D.F90 > CMakeFiles/Gar6more2D.out.dir/lib/bin/Gar6more2D.F90.i

CMakeFiles/Gar6more2D.out.dir/lib/bin/Gar6more2D.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/Gar6more2D.out.dir/lib/bin/Gar6more2D.F90.s"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/mottierr/Codes/Gar6more/Gar6more2d/lib/bin/Gar6more2D.F90 -o CMakeFiles/Gar6more2D.out.dir/lib/bin/Gar6more2D.F90.s

# Object files for target Gar6more2D.out
Gar6more2D_out_OBJECTS = \
"CMakeFiles/Gar6more2D.out.dir/lib/bin/Gar6more2D.F90.o"

# External object files for target Gar6more2D.out
Gar6more2D_out_EXTERNAL_OBJECTS =

Gar6more2D.out: CMakeFiles/Gar6more2D.out.dir/lib/bin/Gar6more2D.F90.o
Gar6more2D.out: CMakeFiles/Gar6more2D.out.dir/build.make
Gar6more2D.out: mod/libmod.a
Gar6more2D.out: lib/libacousacous/libacousacous.a
Gar6more2D.out: lib/libacouselasto/libacouselasto.a
Gar6more2D.out: lib/libacousporo/libacousporo.a
Gar6more2D.out: lib/libelastoelasto/libelastoelasto.a
Gar6more2D.out: lib/libporoporo/libporoporo.a
Gar6more2D.out: lib/libporoelasto/libporoelasto.a
Gar6more2D.out: lib/libgeneral/libgene.a
Gar6more2D.out: /usr/lib/x86_64-linux-gnu/libopenblas.so
Gar6more2D.out: CMakeFiles/Gar6more2D.out.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mottierr/Codes/Gar6more/Gar6more2d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable Gar6more2D.out"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Gar6more2D.out.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Gar6more2D.out.dir/build: Gar6more2D.out
.PHONY : CMakeFiles/Gar6more2D.out.dir/build

CMakeFiles/Gar6more2D.out.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Gar6more2D.out.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Gar6more2D.out.dir/clean

CMakeFiles/Gar6more2D.out.dir/depend:
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mottierr/Codes/Gar6more/Gar6more2d /home/mottierr/Codes/Gar6more/Gar6more2d /home/mottierr/Codes/Gar6more/Gar6more2d/build /home/mottierr/Codes/Gar6more/Gar6more2d/build /home/mottierr/Codes/Gar6more/Gar6more2d/build/CMakeFiles/Gar6more2D.out.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Gar6more2D.out.dir/depend
