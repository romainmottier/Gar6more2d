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
include mod/CMakeFiles/mod.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include mod/CMakeFiles/mod.dir/compiler_depend.make

# Include the progress variables for this target.
include mod/CMakeFiles/mod.dir/progress.make

# Include the compile flags for this target's objects.
include mod/CMakeFiles/mod.dir/flags.make

mod/CMakeFiles/mod.dir/m_const.F90.o: mod/CMakeFiles/mod.dir/flags.make
mod/CMakeFiles/mod.dir/m_const.F90.o: ../mod/m_const.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mottierr/Codes/Gar6more/Gar6more2d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object mod/CMakeFiles/mod.dir/m_const.F90.o"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_const.F90 -o CMakeFiles/mod.dir/m_const.F90.o

mod/CMakeFiles/mod.dir/m_const.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/mod.dir/m_const.F90.i"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_const.F90 > CMakeFiles/mod.dir/m_const.F90.i

mod/CMakeFiles/mod.dir/m_const.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/mod.dir/m_const.F90.s"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_const.F90 -o CMakeFiles/mod.dir/m_const.F90.s

mod/CMakeFiles/mod.dir/m_num.F90.o: mod/CMakeFiles/mod.dir/flags.make
mod/CMakeFiles/mod.dir/m_num.F90.o: ../mod/m_num.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mottierr/Codes/Gar6more/Gar6more2d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object mod/CMakeFiles/mod.dir/m_num.F90.o"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_num.F90 -o CMakeFiles/mod.dir/m_num.F90.o

mod/CMakeFiles/mod.dir/m_num.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/mod.dir/m_num.F90.i"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_num.F90 > CMakeFiles/mod.dir/m_num.F90.i

mod/CMakeFiles/mod.dir/m_num.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/mod.dir/m_num.F90.s"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_num.F90 -o CMakeFiles/mod.dir/m_num.F90.s

mod/CMakeFiles/mod.dir/m_phys.F90.o: mod/CMakeFiles/mod.dir/flags.make
mod/CMakeFiles/mod.dir/m_phys.F90.o: ../mod/m_phys.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mottierr/Codes/Gar6more/Gar6more2d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object mod/CMakeFiles/mod.dir/m_phys.F90.o"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_phys.F90 -o CMakeFiles/mod.dir/m_phys.F90.o

mod/CMakeFiles/mod.dir/m_phys.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/mod.dir/m_phys.F90.i"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_phys.F90 > CMakeFiles/mod.dir/m_phys.F90.i

mod/CMakeFiles/mod.dir/m_phys.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/mod.dir/m_phys.F90.s"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_phys.F90 -o CMakeFiles/mod.dir/m_phys.F90.s

mod/CMakeFiles/mod.dir/m_result.F90.o: mod/CMakeFiles/mod.dir/flags.make
mod/CMakeFiles/mod.dir/m_result.F90.o: ../mod/m_result.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mottierr/Codes/Gar6more/Gar6more2d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building Fortran object mod/CMakeFiles/mod.dir/m_result.F90.o"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_result.F90 -o CMakeFiles/mod.dir/m_result.F90.o

mod/CMakeFiles/mod.dir/m_result.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/mod.dir/m_result.F90.i"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_result.F90 > CMakeFiles/mod.dir/m_result.F90.i

mod/CMakeFiles/mod.dir/m_result.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/mod.dir/m_result.F90.s"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_result.F90 -o CMakeFiles/mod.dir/m_result.F90.s

mod/CMakeFiles/mod.dir/m_sismo.F90.o: mod/CMakeFiles/mod.dir/flags.make
mod/CMakeFiles/mod.dir/m_sismo.F90.o: ../mod/m_sismo.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mottierr/Codes/Gar6more/Gar6more2d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building Fortran object mod/CMakeFiles/mod.dir/m_sismo.F90.o"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_sismo.F90 -o CMakeFiles/mod.dir/m_sismo.F90.o

mod/CMakeFiles/mod.dir/m_sismo.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/mod.dir/m_sismo.F90.i"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_sismo.F90 > CMakeFiles/mod.dir/m_sismo.F90.i

mod/CMakeFiles/mod.dir/m_sismo.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/mod.dir/m_sismo.F90.s"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_sismo.F90 -o CMakeFiles/mod.dir/m_sismo.F90.s

mod/CMakeFiles/mod.dir/m_source.F90.o: mod/CMakeFiles/mod.dir/flags.make
mod/CMakeFiles/mod.dir/m_source.F90.o: ../mod/m_source.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mottierr/Codes/Gar6more/Gar6more2d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building Fortran object mod/CMakeFiles/mod.dir/m_source.F90.o"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_source.F90 -o CMakeFiles/mod.dir/m_source.F90.o

mod/CMakeFiles/mod.dir/m_source.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/mod.dir/m_source.F90.i"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_source.F90 > CMakeFiles/mod.dir/m_source.F90.i

mod/CMakeFiles/mod.dir/m_source.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/mod.dir/m_source.F90.s"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/mottierr/Codes/Gar6more/Gar6more2d/mod/m_source.F90 -o CMakeFiles/mod.dir/m_source.F90.s

# Object files for target mod
mod_OBJECTS = \
"CMakeFiles/mod.dir/m_const.F90.o" \
"CMakeFiles/mod.dir/m_num.F90.o" \
"CMakeFiles/mod.dir/m_phys.F90.o" \
"CMakeFiles/mod.dir/m_result.F90.o" \
"CMakeFiles/mod.dir/m_sismo.F90.o" \
"CMakeFiles/mod.dir/m_source.F90.o"

# External object files for target mod
mod_EXTERNAL_OBJECTS =

mod/libmod.a: mod/CMakeFiles/mod.dir/m_const.F90.o
mod/libmod.a: mod/CMakeFiles/mod.dir/m_num.F90.o
mod/libmod.a: mod/CMakeFiles/mod.dir/m_phys.F90.o
mod/libmod.a: mod/CMakeFiles/mod.dir/m_result.F90.o
mod/libmod.a: mod/CMakeFiles/mod.dir/m_sismo.F90.o
mod/libmod.a: mod/CMakeFiles/mod.dir/m_source.F90.o
mod/libmod.a: mod/CMakeFiles/mod.dir/build.make
mod/libmod.a: mod/CMakeFiles/mod.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mottierr/Codes/Gar6more/Gar6more2d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking Fortran static library libmod.a"
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && $(CMAKE_COMMAND) -P CMakeFiles/mod.dir/cmake_clean_target.cmake
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mod.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
mod/CMakeFiles/mod.dir/build: mod/libmod.a
.PHONY : mod/CMakeFiles/mod.dir/build

mod/CMakeFiles/mod.dir/clean:
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod && $(CMAKE_COMMAND) -P CMakeFiles/mod.dir/cmake_clean.cmake
.PHONY : mod/CMakeFiles/mod.dir/clean

mod/CMakeFiles/mod.dir/depend:
	cd /home/mottierr/Codes/Gar6more/Gar6more2d/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mottierr/Codes/Gar6more/Gar6more2d /home/mottierr/Codes/Gar6more/Gar6more2d/mod /home/mottierr/Codes/Gar6more/Gar6more2d/build /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod /home/mottierr/Codes/Gar6more/Gar6more2d/build/mod/CMakeFiles/mod.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : mod/CMakeFiles/mod.dir/depend

