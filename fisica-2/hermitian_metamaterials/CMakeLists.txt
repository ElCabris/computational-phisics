# =============================================================================
# CMakeLists.txt — hermetic-bands
# Photonic band structure of dispersive metamaterials via Hermitian eigenvalue
# formulation.  Raman & Fan, PRL 104, 087401 (2010).
#
# Build:
#   mkdir build && cd build
#   cmake -DCMAKE_BUILD_TYPE=Release ..
#   make -j
# =============================================================================

cmake_minimum_required(VERSION 3.16)
project(hermetic_bands
    VERSION 1.0.0
    DESCRIPTION "Hermitian photonic band structure solver for dispersive metamaterials"
    LANGUAGES C CXX Fortran
)

# ---------------------------------------------------------------------------
# Language standards
# ---------------------------------------------------------------------------
set(CMAKE_C_STANDARD   17)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD   17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_Fortran_FORMAT FREE)

# ---------------------------------------------------------------------------
# Default build type
# ---------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ---------------------------------------------------------------------------
# Compiler flags
# ---------------------------------------------------------------------------

# --- C++ ---
set(CXX_COMMON_FLAGS "-Wall -Wextra -Wpedantic -march=native -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE  "-O3 ${CXX_COMMON_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG    "-O0 -g3 -DDEBUG ${CXX_COMMON_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g ${CXX_COMMON_FLAGS}")

# --- Fortran ---
set(FC_COMMON_FLAGS "-Wall -Wextra")
set(CMAKE_Fortran_FLAGS_RELEASE  "-O3 ${FC_COMMON_FLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG    "-O0 -g3 -fcheck=all -fbacktrace ${FC_COMMON_FLAGS}")

# ---------------------------------------------------------------------------
# Find required libraries: BLAS, LAPACK, ARPACK
# ---------------------------------------------------------------------------

# BLAS
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
    message(STATUS "BLAS found: ${BLAS_LIBRARIES}")
else()
    message(FATAL_ERROR "BLAS not found. Install libblas-dev or openblas-dev.")
endif()

# LAPACK
find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
    message(STATUS "LAPACK found: ${LAPACK_LIBRARIES}")
else()
    message(FATAL_ERROR "LAPACK not found. Install liblapack-dev.")
endif()

# ARPACK — not a standard CMake module; search manually.
# On most distros: libparpack2-dev / libarpack2-dev.
# We look for the Fortran interface library.
find_library(ARPACK_LIB
    NAMES arpack parpack arpack_ng
    PATHS
        /usr/lib
        /usr/lib/x86_64-linux-gnu
        /usr/lib/aarch64-linux-gnu
        /usr/local/lib
        /usr/local/lib64
        /opt/local/lib
        /opt/homebrew/lib
    PATH_SUFFIXES lib lib64
    DOC "ARPACK library"
)
if(NOT ARPACK_LIB)
    message(FATAL_ERROR
        "ARPACK not found. Install:\n"
        "  Arch Linux    : sudo pacman -S arpack\n"
        "  Debian/Ubuntu : sudo apt install libarpack2-dev\n"
        "  Fedora        : sudo dnf install arpack-devel\n"
        "  openSUSE      : sudo zypper install arpack-devel\n"
        "Then re-run: cmake .."
    )
else()
    message(STATUS "ARPACK found: ${ARPACK_LIB}")
endif()

# Also find gfortran runtime (needed when mixing C++ and Fortran objects)
find_package(OpenMP QUIET)   # optional, for future OMP parallelism

# ---------------------------------------------------------------------------
# Fortran module files go to a dedicated directory
# ---------------------------------------------------------------------------
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_modules)

# ---------------------------------------------------------------------------
# Source lists
# ---------------------------------------------------------------------------

set(FORTRAN_SOURCES
    fortran/precision_kinds.f90
    fortran/arpack_interface.f90
    fortran/lanczos_helpers.f90
)

set(CXX_SOURCES
    src/utils.cpp
    src/sparse.cpp
    src/yee_grid.cpp
    src/operators.cpp
    src/matrices.cpp
    src/eigensolver.cpp
    src/perturbation.cpp
    src/output.cpp
    src/main.cpp
)

# ---------------------------------------------------------------------------
# Build Fortran objects as a static library first.
# This is necessary so that CMake generates the .mod files before the C++
# translation units that include Fortran-callable wrappers are compiled.
# ---------------------------------------------------------------------------
add_library(hermetic_fortran STATIC ${FORTRAN_SOURCES})
target_include_directories(hermetic_fortran PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# ---------------------------------------------------------------------------
# Main executable
# ---------------------------------------------------------------------------
add_executable(hermetic-bands ${CXX_SOURCES})

target_include_directories(hermetic-bands PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Link order matters: our Fortran lib first, then ARPACK (which calls LAPACK/BLAS),
# then LAPACK, then BLAS, then the system math library.
target_link_libraries(hermetic-bands
    hermetic_fortran
    ${ARPACK_LIB}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
    m                       # libm — standard math (sqrt, etc.)
    gfortran                # gfortran runtime for Fortran intrinsics
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(hermetic-bands OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP enabled.")
endif()

# ---------------------------------------------------------------------------
# Install target
# ---------------------------------------------------------------------------
install(TARGETS hermetic-bands DESTINATION bin)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message(STATUS "=== hermetic-bands configuration summary ===")
message(STATUS "  C   standard  : C${CMAKE_C_STANDARD}")
message(STATUS "  C++ standard  : C++${CMAKE_CXX_STANDARD}")
message(STATUS "  Build type    : ${CMAKE_BUILD_TYPE}")
message(STATUS "  BLAS          : ${BLAS_LIBRARIES}")
message(STATUS "  LAPACK        : ${LAPACK_LIBRARIES}")
message(STATUS "  ARPACK        : ${ARPACK_LIB}")
message(STATUS "=============================================")
