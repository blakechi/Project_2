cmake_minimum_required(VERSION 3.1)
project(Project_2 VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)
set(ALL_INCLUDE_DIRECTORIES)

# === add executable
add_executable(project
    "main.cpp"
    "src/Color.hpp"
    "src/Color.cpp"
    "src/ColorMap.hpp"
    "src/DataSource.hpp"
    "src/ScatterPoint.hpp"
    "src/Interpolation.hpp"
    "src/RandNumGenerator.hpp"
    "src/Timer.hpp"
    "src/Timer.cpp"
    "src/utils.hpp"
    "src/utils.cpp"
    "src/Vector2.hpp"
    "src/Vector3.hpp"
    "src/Vector3.cpp"
    "dependence/nanoflann/nanoflann.hpp"
    "dependence/nanoflann/util.h"
    "dependence/reader.hpp"
    "dependence/reader.cpp"
    )
target_include_directories(project
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/dependence/"
    PUBLIC ${ALL_INCLUDE_DIRECTORIES}
    )


# === libraries inside this project
# glm
if(EXISTS ${PROJECT_SOURCE_DIR}/dependence/glm)
    add_subdirectory(${PROJECT_SOURCE_DIR}/dependence/glm)
    target_include_directories(project PUBLIC ${PROJECT_SOURCE_DIR}/dependence/glm)
else()
    message(STATUS "glm doesn't exist")
endif()

# rapidjson
if(EXISTS ${PROJECT_SOURCE_DIR}/dependence/rapidjson)
    add_subdirectory(${PROJECT_SOURCE_DIR}/dependence/rapidjson)
    target_include_directories(project PUBLIC ${PROJECT_SOURCE_DIR}/dependence/rapidjson/include)
else()
    message(STATUS "RapidJSON doesn't exist")
endif()

# === self-defined target
add_custom_target(run
    COMMENT "Running..."
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS project
    )
add_custom_command(TARGET run
    COMMAND ./project
    DEPENDS project
    )


