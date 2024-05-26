cmake_minimum_required(VERSION 3.25)
set(CMAKE_CXX_STANDARD 17)

enable_testing()

project(hmshtest)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Using clang/gcc compiler flags")
    SET(BASE_CXX_FLAGS "-std=c++20 -Wall -Wextra -g3")
    SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-unused-private-field")
    SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        message("Setting clang-specific options")
        SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
        SET(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address -fno-limit-debug-info")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        message("Setting gcc-specific options")
        SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
        SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
        SET(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address")
    endif()

    SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3")
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG(-march=native  COMPILER_SUPPORTS_MARCH_NATIVE)
    if(COMPILER_SUPPORTS_MARCH_NATIVE)
        set(MARCH_NATIVE "-march=native")
    else()
        set(MARCH_NATIVE "")
    endif()
    SET(CMAKE_CXX_FLAGS_RELEASE "${MARCH_NATIVE} -O3 -DNDEBUG")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message("Using Visual Studio compiler flags")
    set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
    set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
    SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

    add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
    add_definitions (-DNOMINMAX)
else()
    message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()


### Download and unpack googletest at configure time
# (from the googletest docs)
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(
        COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(
        COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

add_subdirectory(
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-src
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)

set_property(TARGET gtest PROPERTY CXX_STANDARD 17)

include(FetchContent)
FetchContent_Declare(
        libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG v2.5.0
)
FetchContent_MakeAvailable(libigl)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. "hmsh")
add_executable(${PROJECT_NAME} test.cpp)
target_link_libraries(${PROJECT_NAME} gtest_main igl::core hmsh)

