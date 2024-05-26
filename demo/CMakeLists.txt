cmake_minimum_required(VERSION 3.28)
project(hmshdemo)

set(CMAKE_CXX_STANDARD 17)

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Using clang/gcc compiler flags")
    set(BASE_CXX_FLAGS "-std=c++17 -Wall -Wextra -g3")
    set(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces")
    set(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        message("Setting clang-specific options")
        set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=5 -fcolor-diagnostics")
        set(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address -fno-limit-debug-info")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
        message("Setting gcc-specific options")
        set(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
    endif()

    set(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3")
    set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
    if(COMPILER_SUPPORTS_MARCH_NATIVE)
        set(MARCH_NATIVE "-march=native")
    else()
        set(MARCH_NATIVE "")
    endif()
    set(CMAKE_CXX_FLAGS_RELEASE "${MARCH_NATIVE} -O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
else()
    message(FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()

include(FetchContent)
FetchContent_Declare(
        libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG v2.5.0
)
FetchContent_MakeAvailable(libigl)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../polyscope/ "polyscope")
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../ "hmsh")
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE polyscope igl::core hmsh)
