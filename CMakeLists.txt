cmake_minimum_required(VERSION 3.5)
project(ESGI140)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-O3 -march=native)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/") 
#include(FindOpenMP)
#if(OPENMP_FOUND)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif(OPENMP_FOUND)


find_package(Boost REQUIRED) 

add_library(${PROJECT_NAME} STATIC
    src/particlemodel.cpp
    src/measure.cpp
)

set_target_properties(${PROJECT_NAME}
    PROPERTIES POSITION_INDEPENDENT_CODE ON
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON)

target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS})

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

