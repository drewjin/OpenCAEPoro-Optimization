# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /public4/home/sc56808/zlj/parmetis-4.0.3/metis

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64

# Include any dependencies generated for this target.
include programs/CMakeFiles/mpmetis.dir/depend.make

# Include the progress variables for this target.
include programs/CMakeFiles/mpmetis.dir/progress.make

# Include the compile flags for this target's objects.
include programs/CMakeFiles/mpmetis.dir/flags.make

programs/CMakeFiles/mpmetis.dir/mpmetis.c.o: programs/CMakeFiles/mpmetis.dir/flags.make
programs/CMakeFiles/mpmetis.dir/mpmetis.c.o: ../../programs/mpmetis.c
	$(CMAKE_COMMAND) -E cmake_progress_report /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object programs/CMakeFiles/mpmetis.dir/mpmetis.c.o"
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs && /public4/soft/oneAPI/2022.1/mpi/2021.5.0/bin/mpiicc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/mpmetis.dir/mpmetis.c.o   -c /public4/home/sc56808/zlj/parmetis-4.0.3/metis/programs/mpmetis.c

programs/CMakeFiles/mpmetis.dir/mpmetis.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mpmetis.dir/mpmetis.c.i"
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs && /public4/soft/oneAPI/2022.1/mpi/2021.5.0/bin/mpiicc  $(C_DEFINES) $(C_FLAGS) -E /public4/home/sc56808/zlj/parmetis-4.0.3/metis/programs/mpmetis.c > CMakeFiles/mpmetis.dir/mpmetis.c.i

programs/CMakeFiles/mpmetis.dir/mpmetis.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mpmetis.dir/mpmetis.c.s"
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs && /public4/soft/oneAPI/2022.1/mpi/2021.5.0/bin/mpiicc  $(C_DEFINES) $(C_FLAGS) -S /public4/home/sc56808/zlj/parmetis-4.0.3/metis/programs/mpmetis.c -o CMakeFiles/mpmetis.dir/mpmetis.c.s

programs/CMakeFiles/mpmetis.dir/mpmetis.c.o.requires:
.PHONY : programs/CMakeFiles/mpmetis.dir/mpmetis.c.o.requires

programs/CMakeFiles/mpmetis.dir/mpmetis.c.o.provides: programs/CMakeFiles/mpmetis.dir/mpmetis.c.o.requires
	$(MAKE) -f programs/CMakeFiles/mpmetis.dir/build.make programs/CMakeFiles/mpmetis.dir/mpmetis.c.o.provides.build
.PHONY : programs/CMakeFiles/mpmetis.dir/mpmetis.c.o.provides

programs/CMakeFiles/mpmetis.dir/mpmetis.c.o.provides.build: programs/CMakeFiles/mpmetis.dir/mpmetis.c.o

programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o: programs/CMakeFiles/mpmetis.dir/flags.make
programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o: ../../programs/cmdline_mpmetis.c
	$(CMAKE_COMMAND) -E cmake_progress_report /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o"
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs && /public4/soft/oneAPI/2022.1/mpi/2021.5.0/bin/mpiicc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o   -c /public4/home/sc56808/zlj/parmetis-4.0.3/metis/programs/cmdline_mpmetis.c

programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.i"
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs && /public4/soft/oneAPI/2022.1/mpi/2021.5.0/bin/mpiicc  $(C_DEFINES) $(C_FLAGS) -E /public4/home/sc56808/zlj/parmetis-4.0.3/metis/programs/cmdline_mpmetis.c > CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.i

programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.s"
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs && /public4/soft/oneAPI/2022.1/mpi/2021.5.0/bin/mpiicc  $(C_DEFINES) $(C_FLAGS) -S /public4/home/sc56808/zlj/parmetis-4.0.3/metis/programs/cmdline_mpmetis.c -o CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.s

programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o.requires:
.PHONY : programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o.requires

programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o.provides: programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o.requires
	$(MAKE) -f programs/CMakeFiles/mpmetis.dir/build.make programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o.provides.build
.PHONY : programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o.provides

programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o.provides.build: programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o

programs/CMakeFiles/mpmetis.dir/io.c.o: programs/CMakeFiles/mpmetis.dir/flags.make
programs/CMakeFiles/mpmetis.dir/io.c.o: ../../programs/io.c
	$(CMAKE_COMMAND) -E cmake_progress_report /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object programs/CMakeFiles/mpmetis.dir/io.c.o"
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs && /public4/soft/oneAPI/2022.1/mpi/2021.5.0/bin/mpiicc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/mpmetis.dir/io.c.o   -c /public4/home/sc56808/zlj/parmetis-4.0.3/metis/programs/io.c

programs/CMakeFiles/mpmetis.dir/io.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mpmetis.dir/io.c.i"
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs && /public4/soft/oneAPI/2022.1/mpi/2021.5.0/bin/mpiicc  $(C_DEFINES) $(C_FLAGS) -E /public4/home/sc56808/zlj/parmetis-4.0.3/metis/programs/io.c > CMakeFiles/mpmetis.dir/io.c.i

programs/CMakeFiles/mpmetis.dir/io.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mpmetis.dir/io.c.s"
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs && /public4/soft/oneAPI/2022.1/mpi/2021.5.0/bin/mpiicc  $(C_DEFINES) $(C_FLAGS) -S /public4/home/sc56808/zlj/parmetis-4.0.3/metis/programs/io.c -o CMakeFiles/mpmetis.dir/io.c.s

programs/CMakeFiles/mpmetis.dir/io.c.o.requires:
.PHONY : programs/CMakeFiles/mpmetis.dir/io.c.o.requires

programs/CMakeFiles/mpmetis.dir/io.c.o.provides: programs/CMakeFiles/mpmetis.dir/io.c.o.requires
	$(MAKE) -f programs/CMakeFiles/mpmetis.dir/build.make programs/CMakeFiles/mpmetis.dir/io.c.o.provides.build
.PHONY : programs/CMakeFiles/mpmetis.dir/io.c.o.provides

programs/CMakeFiles/mpmetis.dir/io.c.o.provides.build: programs/CMakeFiles/mpmetis.dir/io.c.o

programs/CMakeFiles/mpmetis.dir/stat.c.o: programs/CMakeFiles/mpmetis.dir/flags.make
programs/CMakeFiles/mpmetis.dir/stat.c.o: ../../programs/stat.c
	$(CMAKE_COMMAND) -E cmake_progress_report /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object programs/CMakeFiles/mpmetis.dir/stat.c.o"
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs && /public4/soft/oneAPI/2022.1/mpi/2021.5.0/bin/mpiicc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/mpmetis.dir/stat.c.o   -c /public4/home/sc56808/zlj/parmetis-4.0.3/metis/programs/stat.c

programs/CMakeFiles/mpmetis.dir/stat.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mpmetis.dir/stat.c.i"
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs && /public4/soft/oneAPI/2022.1/mpi/2021.5.0/bin/mpiicc  $(C_DEFINES) $(C_FLAGS) -E /public4/home/sc56808/zlj/parmetis-4.0.3/metis/programs/stat.c > CMakeFiles/mpmetis.dir/stat.c.i

programs/CMakeFiles/mpmetis.dir/stat.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mpmetis.dir/stat.c.s"
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs && /public4/soft/oneAPI/2022.1/mpi/2021.5.0/bin/mpiicc  $(C_DEFINES) $(C_FLAGS) -S /public4/home/sc56808/zlj/parmetis-4.0.3/metis/programs/stat.c -o CMakeFiles/mpmetis.dir/stat.c.s

programs/CMakeFiles/mpmetis.dir/stat.c.o.requires:
.PHONY : programs/CMakeFiles/mpmetis.dir/stat.c.o.requires

programs/CMakeFiles/mpmetis.dir/stat.c.o.provides: programs/CMakeFiles/mpmetis.dir/stat.c.o.requires
	$(MAKE) -f programs/CMakeFiles/mpmetis.dir/build.make programs/CMakeFiles/mpmetis.dir/stat.c.o.provides.build
.PHONY : programs/CMakeFiles/mpmetis.dir/stat.c.o.provides

programs/CMakeFiles/mpmetis.dir/stat.c.o.provides.build: programs/CMakeFiles/mpmetis.dir/stat.c.o

# Object files for target mpmetis
mpmetis_OBJECTS = \
"CMakeFiles/mpmetis.dir/mpmetis.c.o" \
"CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o" \
"CMakeFiles/mpmetis.dir/io.c.o" \
"CMakeFiles/mpmetis.dir/stat.c.o"

# External object files for target mpmetis
mpmetis_EXTERNAL_OBJECTS =

programs/mpmetis: programs/CMakeFiles/mpmetis.dir/mpmetis.c.o
programs/mpmetis: programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o
programs/mpmetis: programs/CMakeFiles/mpmetis.dir/io.c.o
programs/mpmetis: programs/CMakeFiles/mpmetis.dir/stat.c.o
programs/mpmetis: programs/CMakeFiles/mpmetis.dir/build.make
programs/mpmetis: libmetis/libmetis.a
programs/mpmetis: programs/CMakeFiles/mpmetis.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable mpmetis"
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mpmetis.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
programs/CMakeFiles/mpmetis.dir/build: programs/mpmetis
.PHONY : programs/CMakeFiles/mpmetis.dir/build

programs/CMakeFiles/mpmetis.dir/requires: programs/CMakeFiles/mpmetis.dir/mpmetis.c.o.requires
programs/CMakeFiles/mpmetis.dir/requires: programs/CMakeFiles/mpmetis.dir/cmdline_mpmetis.c.o.requires
programs/CMakeFiles/mpmetis.dir/requires: programs/CMakeFiles/mpmetis.dir/io.c.o.requires
programs/CMakeFiles/mpmetis.dir/requires: programs/CMakeFiles/mpmetis.dir/stat.c.o.requires
.PHONY : programs/CMakeFiles/mpmetis.dir/requires

programs/CMakeFiles/mpmetis.dir/clean:
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs && $(CMAKE_COMMAND) -P CMakeFiles/mpmetis.dir/cmake_clean.cmake
.PHONY : programs/CMakeFiles/mpmetis.dir/clean

programs/CMakeFiles/mpmetis.dir/depend:
	cd /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /public4/home/sc56808/zlj/parmetis-4.0.3/metis /public4/home/sc56808/zlj/parmetis-4.0.3/metis/programs /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64 /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs /public4/home/sc56808/zlj/parmetis-4.0.3/metis/build/Linux-x86_64/programs/CMakeFiles/mpmetis.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : programs/CMakeFiles/mpmetis.dir/depend

