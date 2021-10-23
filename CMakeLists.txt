cmake_minimum_required(VERSION 3.10)

#project name
project(GridProjector)

option(CXX "enable C++ compilation" ON)
enable_language(CXX)

#specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")
SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")


message("BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

add_executable(${PROJECT_NAME} main.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC
	./
	external/cinolib/external/eigen/
	external/cinolib/include/
	${GUROBI_INCLUDE_DIRS}
)