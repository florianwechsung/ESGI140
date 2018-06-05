cmake_minimum_required(VERSION 3.5)
project(ESGI140)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-O3 -march=native)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/") 

add_definitions("-DEIGEN_DONT_PARALLELIZE")


find_package(Boost REQUIRED) 
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_library(${PROJECT_NAME} STATIC
    src/particlemodel.cpp
)

set_target_properties(${PROJECT_NAME}
    PROPERTIES POSITION_INDEPENDENT_CODE ON
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON)

target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} Eigen3::Eigen)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include ${EIGEN3_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} )

set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)

add_executable(driver_particlemodel driver/driver_particlemodel.cpp)
target_link_libraries(driver_particlemodel ${PROJECT_NAME})
target_include_directories(driver_particlemodel PRIVATE ${PROJECT_SOURCE_DIR}/include)
set_target_properties(driver_particlemodel PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED ON)



add_subdirectory(pybind11)
pybind11_add_module(pyesgi140 pythonbindings/main.cpp)
set_target_properties(pyesgi140 PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED ON)
target_link_libraries(pyesgi140 PRIVATE ${PROJECT_NAME})
target_include_directories(pyesgi140 PRIVATE ${${PROJECT_NAME}_INCLUDE_DIRS})

