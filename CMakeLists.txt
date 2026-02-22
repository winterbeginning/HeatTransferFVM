cmake_minimum_required(VERSION 3.10)
project(HeatTransferFVM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(include)

file(GLOB SOURCES "src/*.cpp")

add_executable(HeatTransferFVM ${SOURCES})
