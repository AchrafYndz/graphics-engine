# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /opt/clion-2021.2.3/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /opt/clion-2021.2.3/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/achraf/CLionProjects/graphics-engine

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/achraf/CLionProjects/graphics-engine/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/engine.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/engine.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/engine.dir/flags.make

CMakeFiles/engine.dir/util/Easy_image/easy_image.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/util/Easy_image/easy_image.cc.o: ../util/Easy_image/easy_image.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/achraf/CLionProjects/graphics-engine/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/engine.dir/util/Easy_image/easy_image.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/util/Easy_image/easy_image.cc.o -c /home/achraf/CLionProjects/graphics-engine/util/Easy_image/easy_image.cc

CMakeFiles/engine.dir/util/Easy_image/easy_image.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/util/Easy_image/easy_image.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/achraf/CLionProjects/graphics-engine/util/Easy_image/easy_image.cc > CMakeFiles/engine.dir/util/Easy_image/easy_image.cc.i

CMakeFiles/engine.dir/util/Easy_image/easy_image.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/util/Easy_image/easy_image.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/achraf/CLionProjects/graphics-engine/util/Easy_image/easy_image.cc -o CMakeFiles/engine.dir/util/Easy_image/easy_image.cc.s

CMakeFiles/engine.dir/engine.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/engine.cc.o: ../engine.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/achraf/CLionProjects/graphics-engine/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/engine.dir/engine.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/engine.cc.o -c /home/achraf/CLionProjects/graphics-engine/engine.cc

CMakeFiles/engine.dir/engine.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/engine.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/achraf/CLionProjects/graphics-engine/engine.cc > CMakeFiles/engine.dir/engine.cc.i

CMakeFiles/engine.dir/engine.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/engine.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/achraf/CLionProjects/graphics-engine/engine.cc -o CMakeFiles/engine.dir/engine.cc.s

CMakeFiles/engine.dir/util/Ini_config/ini_configuration.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/util/Ini_config/ini_configuration.cc.o: ../util/Ini_config/ini_configuration.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/achraf/CLionProjects/graphics-engine/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/engine.dir/util/Ini_config/ini_configuration.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/util/Ini_config/ini_configuration.cc.o -c /home/achraf/CLionProjects/graphics-engine/util/Ini_config/ini_configuration.cc

CMakeFiles/engine.dir/util/Ini_config/ini_configuration.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/util/Ini_config/ini_configuration.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/achraf/CLionProjects/graphics-engine/util/Ini_config/ini_configuration.cc > CMakeFiles/engine.dir/util/Ini_config/ini_configuration.cc.i

CMakeFiles/engine.dir/util/Ini_config/ini_configuration.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/util/Ini_config/ini_configuration.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/achraf/CLionProjects/graphics-engine/util/Ini_config/ini_configuration.cc -o CMakeFiles/engine.dir/util/Ini_config/ini_configuration.cc.s

CMakeFiles/engine.dir/util/Color.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/util/Color.cc.o: ../util/Color.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/achraf/CLionProjects/graphics-engine/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/engine.dir/util/Color.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/util/Color.cc.o -c /home/achraf/CLionProjects/graphics-engine/util/Color.cc

CMakeFiles/engine.dir/util/Color.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/util/Color.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/achraf/CLionProjects/graphics-engine/util/Color.cc > CMakeFiles/engine.dir/util/Color.cc.i

CMakeFiles/engine.dir/util/Color.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/util/Color.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/achraf/CLionProjects/graphics-engine/util/Color.cc -o CMakeFiles/engine.dir/util/Color.cc.s

CMakeFiles/engine.dir/util/Line2D.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/util/Line2D.cc.o: ../util/Line2D.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/achraf/CLionProjects/graphics-engine/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/engine.dir/util/Line2D.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/util/Line2D.cc.o -c /home/achraf/CLionProjects/graphics-engine/util/Line2D.cc

CMakeFiles/engine.dir/util/Line2D.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/util/Line2D.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/achraf/CLionProjects/graphics-engine/util/Line2D.cc > CMakeFiles/engine.dir/util/Line2D.cc.i

CMakeFiles/engine.dir/util/Line2D.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/util/Line2D.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/achraf/CLionProjects/graphics-engine/util/Line2D.cc -o CMakeFiles/engine.dir/util/Line2D.cc.s

CMakeFiles/engine.dir/util/Point2D.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/util/Point2D.cc.o: ../util/Point2D.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/achraf/CLionProjects/graphics-engine/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/engine.dir/util/Point2D.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/util/Point2D.cc.o -c /home/achraf/CLionProjects/graphics-engine/util/Point2D.cc

CMakeFiles/engine.dir/util/Point2D.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/util/Point2D.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/achraf/CLionProjects/graphics-engine/util/Point2D.cc > CMakeFiles/engine.dir/util/Point2D.cc.i

CMakeFiles/engine.dir/util/Point2D.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/util/Point2D.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/achraf/CLionProjects/graphics-engine/util/Point2D.cc -o CMakeFiles/engine.dir/util/Point2D.cc.s

CMakeFiles/engine.dir/util/Parser/l_parser.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/util/Parser/l_parser.cc.o: ../util/Parser/l_parser.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/achraf/CLionProjects/graphics-engine/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/engine.dir/util/Parser/l_parser.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/util/Parser/l_parser.cc.o -c /home/achraf/CLionProjects/graphics-engine/util/Parser/l_parser.cc

CMakeFiles/engine.dir/util/Parser/l_parser.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/util/Parser/l_parser.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/achraf/CLionProjects/graphics-engine/util/Parser/l_parser.cc > CMakeFiles/engine.dir/util/Parser/l_parser.cc.i

CMakeFiles/engine.dir/util/Parser/l_parser.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/util/Parser/l_parser.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/achraf/CLionProjects/graphics-engine/util/Parser/l_parser.cc -o CMakeFiles/engine.dir/util/Parser/l_parser.cc.s

CMakeFiles/engine.dir/util/Face.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/util/Face.cc.o: ../util/Face.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/achraf/CLionProjects/graphics-engine/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/engine.dir/util/Face.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/util/Face.cc.o -c /home/achraf/CLionProjects/graphics-engine/util/Face.cc

CMakeFiles/engine.dir/util/Face.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/util/Face.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/achraf/CLionProjects/graphics-engine/util/Face.cc > CMakeFiles/engine.dir/util/Face.cc.i

CMakeFiles/engine.dir/util/Face.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/util/Face.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/achraf/CLionProjects/graphics-engine/util/Face.cc -o CMakeFiles/engine.dir/util/Face.cc.s

CMakeFiles/engine.dir/util/Figure.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/util/Figure.cc.o: ../util/Figure.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/achraf/CLionProjects/graphics-engine/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/engine.dir/util/Figure.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/util/Figure.cc.o -c /home/achraf/CLionProjects/graphics-engine/util/Figure.cc

CMakeFiles/engine.dir/util/Figure.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/util/Figure.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/achraf/CLionProjects/graphics-engine/util/Figure.cc > CMakeFiles/engine.dir/util/Figure.cc.i

CMakeFiles/engine.dir/util/Figure.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/util/Figure.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/achraf/CLionProjects/graphics-engine/util/Figure.cc -o CMakeFiles/engine.dir/util/Figure.cc.s

CMakeFiles/engine.dir/util/Vector3D/vector3d.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/util/Vector3D/vector3d.cc.o: ../util/Vector3D/vector3d.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/achraf/CLionProjects/graphics-engine/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/engine.dir/util/Vector3D/vector3d.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/util/Vector3D/vector3d.cc.o -c /home/achraf/CLionProjects/graphics-engine/util/Vector3D/vector3d.cc

CMakeFiles/engine.dir/util/Vector3D/vector3d.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/util/Vector3D/vector3d.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/achraf/CLionProjects/graphics-engine/util/Vector3D/vector3d.cc > CMakeFiles/engine.dir/util/Vector3D/vector3d.cc.i

CMakeFiles/engine.dir/util/Vector3D/vector3d.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/util/Vector3D/vector3d.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/achraf/CLionProjects/graphics-engine/util/Vector3D/vector3d.cc -o CMakeFiles/engine.dir/util/Vector3D/vector3d.cc.s

# Object files for target engine
engine_OBJECTS = \
"CMakeFiles/engine.dir/util/Easy_image/easy_image.cc.o" \
"CMakeFiles/engine.dir/engine.cc.o" \
"CMakeFiles/engine.dir/util/Ini_config/ini_configuration.cc.o" \
"CMakeFiles/engine.dir/util/Color.cc.o" \
"CMakeFiles/engine.dir/util/Line2D.cc.o" \
"CMakeFiles/engine.dir/util/Point2D.cc.o" \
"CMakeFiles/engine.dir/util/Parser/l_parser.cc.o" \
"CMakeFiles/engine.dir/util/Face.cc.o" \
"CMakeFiles/engine.dir/util/Figure.cc.o" \
"CMakeFiles/engine.dir/util/Vector3D/vector3d.cc.o"

# External object files for target engine
engine_EXTERNAL_OBJECTS =

engine: CMakeFiles/engine.dir/util/Easy_image/easy_image.cc.o
engine: CMakeFiles/engine.dir/engine.cc.o
engine: CMakeFiles/engine.dir/util/Ini_config/ini_configuration.cc.o
engine: CMakeFiles/engine.dir/util/Color.cc.o
engine: CMakeFiles/engine.dir/util/Line2D.cc.o
engine: CMakeFiles/engine.dir/util/Point2D.cc.o
engine: CMakeFiles/engine.dir/util/Parser/l_parser.cc.o
engine: CMakeFiles/engine.dir/util/Face.cc.o
engine: CMakeFiles/engine.dir/util/Figure.cc.o
engine: CMakeFiles/engine.dir/util/Vector3D/vector3d.cc.o
engine: CMakeFiles/engine.dir/build.make
engine: CMakeFiles/engine.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/achraf/CLionProjects/graphics-engine/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Linking CXX executable engine"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/engine.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/engine.dir/build: engine
.PHONY : CMakeFiles/engine.dir/build

CMakeFiles/engine.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/engine.dir/cmake_clean.cmake
.PHONY : CMakeFiles/engine.dir/clean

CMakeFiles/engine.dir/depend:
	cd /home/achraf/CLionProjects/graphics-engine/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/achraf/CLionProjects/graphics-engine /home/achraf/CLionProjects/graphics-engine /home/achraf/CLionProjects/graphics-engine/cmake-build-debug /home/achraf/CLionProjects/graphics-engine/cmake-build-debug /home/achraf/CLionProjects/graphics-engine/cmake-build-debug/CMakeFiles/engine.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/engine.dir/depend

