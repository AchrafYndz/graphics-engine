cmake_minimum_required(VERSION 3.5)

project(engine)
############################################################
# Set compiler flags
############################################################
# If supported by your compiler, you can add the -Wall, -Wextra, –fstack-protector-all and -g3 flags here.
set(OWN_GXX_FLAGS "-std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OWN_GXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OWN_GXX_FLAGS}")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OWN_GXX_FLAGS}")

############################################################
# List all sources
############################################################
set(engine_sources
		util/Easy_image/easy_image.cc
		engine.cc
		util/Ini_config/ini_configuration.cc
		util/Color.cc
		util/Line2D.cc
		util/Point2D.cc
		util/Parser/l_parser.cc
		util/Face.cc
		util/Face.h
		util/Figure.cc
		util/Figure.h
		util/Vector3D/vector3d.h
		util/Vector3D/vector3d.cc)

############################################################
# Create an executable
############################################################
set(exe_name "engine")
add_executable( ${exe_name} ${engine_sources} )
install( TARGETS ${exe_name} DESTINATION ${PROJECT_SOURCE_DIR}/ )
