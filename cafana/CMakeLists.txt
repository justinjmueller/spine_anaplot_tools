cmake_minimum_required(VERSION 3.14)
project(example LANGUAGES CXX)

# Set C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set optimization flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
#set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Tree)
include(${ROOT_USE_FILE})

# Environment variables and paths
set(SRPROXY_INC "$ENV{SRPROXY_INC}")
set(OSCLIB_INC "$ENV{OSCLIB_INC}")
set(OSCLIB_LIB "$ENV{OSCLIB_LIB}")
set(EIGEN_INC "$ENV{EIGEN_INC}")
set(SBNANA_INC "$ENV{MRB_INSTALL}/sbnana/$ENV{SBNANA_VERSION}/include")
set(SBNANAOBJ_INC "$ENV{MRB_INSTALL}/sbnanaobj/$ENV{SBNANAOBJ_VERSION}/include")
set(SBNANA_LIB "$ENV{MRB_INSTALL}/sbnana/$ENV{SBNANA_VERSION}/slf7.x86_64.e26.prof/lib")
set(SBNANAOBJ_LIB "$ENV{MRB_INSTALL}/sbnanaobj/$ENV{SBNANAOBJ_VERSION}/slf7.x86_64.e26.prof/lib")

# Add a common library
add_library(common INTERFACE)

# Include directories
target_include_directories(common INTERFACE
    ${ROOT_INCLUDE_DIRS}
    ${SRPROXY_INC}
    ${OSCLIB_INC}
    ${EIGEN_INC}
    ${SBNANA_INC}
    ${SBNANAOBJ_INC}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link directories
target_link_directories(common INTERFACE
    ${OSCLIB_LIB}
    ${SBNANA_LIB}
    ${SBNANAOBJ_LIB}
)

# Link required libraries
target_link_libraries(common INTERFACE
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Tree
    Minuit2
    sbnanaobj_StandardRecordProxy
    OscLib
    CAFAnaCore
    CAFAnaUnfold
    CAFAnaVars
    CAFAnaCuts
    CAFAnaExtrap
    CAFAnaExperiment
    CAFAnaAnalysis
)

target_compile_features(common INTERFACE cxx_std_17)

add_executable(example example.cc)
target_link_libraries(example PRIVATE common)

add_executable(muon2024 muon2024.cc)
target_link_libraries(muon2024 PRIVATE common)

add_executable(spine_benchmark spine_benchmark.cc)
target_link_libraries(spine_benchmark PRIVATE common)

