cmake_minimum_required(VERSION 3.28)
project(VaspUtils VERSION 0.1 LANGUAGES C CXX)

include(FetchContent)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON) 

# Put all executables into bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# ===== spglib (vendored, static) =====
FetchContent_Declare(
    spglib
    GIT_REPOSITORY https://github.com/spglib/spglib.git
    GIT_TAG v2.3.0
)

set(SPGLIB_WITH_Fortran OFF CACHE BOOL "" FORCE)
set(SPGLIB_WITH_Python  OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(spglib)

#file(GLOB SPGLIB_SOURCES "external/spglib/src/*.c")

#add_library(spglib STATIC ${SPGLIB_SOURCES})


#add_library(spglib STATIC
#    external/spglib/src/spglib.c
#    external/spglib/src/cell.c
#    external/spglib/src/debug.c
#    external/spglib/src/delaunay.c
#    external/spglib/src/hall_symbol.c
#    external/spglib/src/kgrid.c
#    external/spglib/src/kpoint.c
#    external/spglib/src/mathfunc.c
#    external/spglib/src/niggli.c
#    external/spglib/src/pointgroup.c
#    external/spglib/src/primitive.c
#    external/spglib/src/refinement.c
#    external/spglib/src/site_symmetry.c
#    external/spglib/src/spacegroup.c
#    external/spglib/src/spin.c
#    external/spglib/src/symmetry.c
#    external/spglib/src/arithmetic.c
#    external/spglib/src/determination.c
#    external/spglib/src/magnetic_spacegroup.c
#    external/spglib/src/msg_database.c
#    external/spglib/src/overlap.c
#    external/spglib/src/sitesym_database.c
#    external/spglib/src/spg_database.c    
#)

#target_include_directories(spglib PUBLIC
#    external/spglib/include
#)


# ===== Core library (no spglib) =====
add_library(vasp_core
    src/poscar_file.cpp
    src/random_utility.cpp
    src/matrix_operations.cpp
)

# Headers exported by the core library
target_include_directories(vasp_core PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

# ===== spglib-based helpers =====
add_library(vasp_spglib
    src/symmetry.cpp
)

target_include_directories(vasp_spglib PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(vasp_spglib PUBLIC symspg vasp_core)

# ===== Standard utilities (no spglib) =====
add_executable(poscar_atom_displace src/poscar_atom_displacement.cpp)
target_link_libraries(poscar_atom_displace PRIVATE vasp_core)

add_executable(poscar_d2c src/poscar_d2c.cpp)
target_link_libraries(poscar_d2c PRIVATE vasp_core)

add_executable(poscar_c2d src/poscar_c2d.cpp)
target_link_libraries(poscar_c2d PRIVATE vasp_core)

# ===== spglib utility =====
add_executable(poscar_symmetry src/poscar_symmetry.cpp)
target_link_libraries(poscar_symmetry PRIVATE vasp_spglib)

