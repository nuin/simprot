cmake_minimum_required(VERSION 3.20)
project(simprot VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Werror)
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files
set(CORE_SOURCES
    src/core/random.cpp
    src/core/config.cpp
)

set(TREE_SOURCES
    src/tree/tree_node.cpp
    src/tree/tree_parser.cpp
)

set(SEQUENCE_SOURCES
    src/sequence/sequence_list.cpp
    src/sequence/alignment.cpp
)

set(EVOLUTION_SOURCES
    src/evolution/substitution_matrix.cpp
    src/evolution/indel_model.cpp
    src/evolution/evolver.cpp
)

set(IO_SOURCES
    src/io/newick_reader.cpp
    src/io/fasta_writer.cpp
)

# Library
add_library(simprot_lib STATIC
    ${CORE_SOURCES}
    ${TREE_SOURCES}
    ${SEQUENCE_SOURCES}
    ${EVOLUTION_SOURCES}
    ${IO_SOURCES}
)

target_include_directories(simprot_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

# Main executable
add_executable(simprot src/main.cpp)
target_link_libraries(simprot PRIVATE simprot_lib)

# Testing with Catch2 (optional)
option(BUILD_TESTS "Build unit tests" ON)
if(BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.4.0
    )
    FetchContent_MakeAvailable(Catch2)

    add_executable(simprot_tests
        tests/test_main.cpp
        tests/test_random.cpp
    )
    target_link_libraries(simprot_tests PRIVATE simprot_lib Catch2::Catch2WithMain)

    include(CTest)
    include(Catch)
    catch_discover_tests(simprot_tests)
endif()

# Install
install(TARGETS simprot DESTINATION bin)
install(DIRECTORY include/simprot DESTINATION include)
