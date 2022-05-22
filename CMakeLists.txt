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
		src/engine.cc
		util/Ini_config/ini_configuration.cc
		src/lib/Color.cc
		src/lib/Line2D.cc
		src/lib/Point2D.cc
		util/Parser/l_parser.cc
		src/lib/Face.cc
		src/lib/Face.h
		src/lib/Figure.cc
		src/lib/Figure.h
		util/Vector3D/vector3d.h
		util/Vector3D/vector3d.cc
		src/include/LSystems2D.cc
		src/include/LSystems2D.h
		src/include/LSystems3D.cc
		src/include/LSystems3D.h
		src/include/Figures3D.cc
		src/include/Figures3D.h
		src/include/Intro.cc
		src/include/Intro.h
		src/include/ZBufferWireframes.cc
		src/include/ZBufferWireframes.h
		src/include/Logic3D.cc
		src/include/Logic3D.h
		src/include/ZBufferTriangles.cc
		src/include/ZBufferTriangles.h
		src/include/Fractals3D.cc
		src/include/Fractals3D.h
        src/lib/Light.cc src/lib/Light.h src/lib/InfLight.cc src/lib/InfLight.h src/lib/PointLight.cc src/lib/PointLight.h)

############################################################
# Create an executable
############################################################
set(exe_name "engine")
add_executable( ${exe_name} ${engine_sources} )
install( TARGETS ${exe_name} DESTINATION ${PROJECT_SOURCE_DIR}/ )
