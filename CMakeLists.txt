cmake_minimum_required(VERSION 3.22)
project(VaspUtils VERSION 0.2.2 LANGUAGES C CXX)

include(FetchContent)

# Required so spglib's static lib can be linked into shared libraries
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ===== CLI11 =====
FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.4.2
)
FetchContent_MakeAvailable(cli11)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# Put all executables into bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

#===== gemmi (header-only, CIF / crystallography) =====
FetchContent_Declare(
    gemmi
    GIT_REPOSITORY https://github.com/project-gemmi/gemmi.git
    GIT_TAG v0.6.7
)
set(GEMMI_BUILD_GEMMI_PROG OFF CACHE BOOL "" FORCE)
set(GEMMI_INSTALL         OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(gemmi)

# ===== spglib (vendored, static) =====
FetchContent_Declare(
    spglib
    GIT_REPOSITORY https://github.com/spglib/spglib.git
    GIT_TAG v2.3.0
)

set(SPGLIB_WITH_Fortran OFF CACHE BOOL "" FORCE)
set(SPGLIB_WITH_Python  OFF CACHE BOOL "" FORCE)
set(SPGLIB_WITH_TESTS   OFF CACHE BOOL "" FORCE)

set(BUILD_TESTING OFF)
FetchContent_MakeAvailable(spglib)
set(BUILD_TESTING ON)

# ===== BLAS / LAPACK =====
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# LAPACKE is usually not exposed as a CMake target,
# so we locate it manually
find_library(LAPACKE_LIB lapacke REQUIRED)


# ===== Core library (no spglib) =====
add_library(vasp_core
    src/poscar_file.cpp
    src/random_utility.cpp
)

# Headers exported by the core library
target_include_directories(vasp_core PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

target_link_libraries(vasp_core
    PUBLIC
    ${LAPACKE_LIB}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

add_library(eigenval2bands src/eigenval2bands.cpp)
target_include_directories(eigenval2bands PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

# ===== spglib-based helpers =====
add_library(vasp_spglib
    src/symmetry.cpp
)

target_include_directories(vasp_spglib PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${spglib_SOURCE_DIR}/src
)

target_link_libraries(vasp_spglib PUBLIC Spglib::symspg vasp_core)

# ===== Elastic constants library =====
add_library(vasp_elastic_lib src/elastic.cpp)
target_include_directories(vasp_elastic_lib PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)
target_link_libraries(vasp_elastic_lib PUBLIC vasp_core)

# ===== Surface slab library =====
add_library(vasp_surface_lib src/surface.cpp)
target_include_directories(vasp_surface_lib PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)
target_link_libraries(vasp_surface_lib PUBLIC vasp_core)

# ===== Standard utilities (no spglib) =====
add_executable(poscar_atom_displace src/poscar_atom_displacement.cpp)
target_link_libraries(poscar_atom_displace PRIVATE vasp_core CLI11::CLI11)

add_executable(poscar_d2c src/poscar_d2c.cpp)
target_link_libraries(poscar_d2c PRIVATE vasp_core CLI11::CLI11)

add_executable(poscar_c2d src/poscar_c2d.cpp)
target_link_libraries(poscar_c2d PRIVATE vasp_core CLI11::CLI11)

add_executable(poscar2ctrls src/poscar_2ctrls.cpp)
target_link_libraries(poscar2ctrls PRIVATE vasp_core CLI11::CLI11)

add_executable(poscar_supercell src/poscar_supercell.cpp)
target_link_libraries(poscar_supercell PRIVATE vasp_core CLI11::CLI11)

add_executable(poscar_convergence_test src/poscar_convergence_test.cpp)
target_link_libraries(poscar_convergence_test PRIVATE vasp_core CLI11::CLI11)

# ===== k-path helper library =====
add_library(vasp_kpath src/kpath.cpp)
target_include_directories(vasp_kpath PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(vasp_kpath PUBLIC vasp_spglib)

# ===== spglib utility =====
add_executable(poscar_symmetry src/poscar_symmetry.cpp)
target_link_libraries(poscar_symmetry PRIVATE vasp_spglib CLI11::CLI11)

add_executable(poscar_kpath src/poscar_kpath.cpp)
target_link_libraries(poscar_kpath PRIVATE vasp_kpath CLI11::CLI11)

add_executable(poscar2primitive src/poscar_2primitive.cpp)
target_link_libraries(poscar2primitive PRIVATE vasp_spglib CLI11::CLI11)

add_executable(poscar2conventional src/poscar_2conventional.cpp)
target_link_libraries(poscar2conventional PRIVATE vasp_spglib CLI11::CLI11)

# ===== surface slab utility =====
add_executable(poscar_surface src/poscar_surface.cpp)
target_link_libraries(poscar_surface PRIVATE vasp_surface_lib vasp_spglib CLI11::CLI11)

# ===== elastic constants utility =====
add_executable(poscar_elastic_deform src/poscar_elastic_deform.cpp)
target_link_libraries(poscar_elastic_deform PRIVATE vasp_elastic_lib vasp_spglib CLI11::CLI11)

# ===== CIF library =====
add_library(cif_lib src/cif_file.cpp)
target_include_directories(cif_lib PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)
target_link_libraries(cif_lib PUBLIC vasp_core gemmi::headers)

# ===== CIF utilities =====
add_executable(poscar2cif src/poscar_2cif.cpp)
target_link_libraries(poscar2cif PRIVATE cif_lib CLI11::CLI11)

add_executable(cif2poscar src/cif_2poscar.cpp)
target_link_libraries(cif2poscar PRIVATE cif_lib CLI11::CLI11)

# ===== DOSCAR library =====
add_library(doscar_lib src/doscar.cpp)
target_include_directories(doscar_lib PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

# ===== Postprocessing utilities (no spglib) =====
add_executable(vasp_eigenval2bands src/vasp_eigenval2bands.cpp)
target_link_libraries(vasp_eigenval2bands PRIVATE eigenval2bands CLI11::CLI11)

add_executable(vasp_fatbands src/vasp_fatbands.cpp)
target_link_libraries(vasp_fatbands PRIVATE vasp_core eigenval2bands CLI11::CLI11)

add_executable(vasp_doscar2dos src/vasp_doscar2dos.cpp)
target_link_libraries(vasp_doscar2dos PRIVATE doscar_lib CLI11::CLI11)

add_executable(vasp_elastic_fit src/vasp_elastic_fit.cpp)
target_link_libraries(vasp_elastic_fit PRIVATE vasp_elastic_lib CLI11::CLI11)

add_executable(vasp_thermo src/vasp_thermo.cpp)
target_link_libraries(vasp_thermo PRIVATE vasp_elastic_lib CLI11::CLI11)

add_executable(vasp_magnetic src/vasp_magnetic.cpp)
target_link_libraries(vasp_magnetic PRIVATE vasp_core CLI11::CLI11)

add_executable(vasp_convergence_test src/vasp_convergence_test.cpp)
target_link_libraries(vasp_convergence_test PRIVATE vasp_core CLI11::CLI11)


# ===== Tests (GoogleTest) =====
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()
include(GoogleTest)

add_executable(vasp_tests
    tests/test_poscar_io.cpp
    tests/test_coordinate_conversion.cpp
    tests/test_displacement.cpp
)

target_link_libraries(vasp_tests PRIVATE vasp_core GTest::gtest_main)

# Make test data available relative to the test binary
target_compile_definitions(vasp_tests PRIVATE
    TEST_DATA_DIR="${PROJECT_SOURCE_DIR}"
)

gtest_discover_tests(vasp_tests)

add_executable(vasp_kpath_tests
    tests/test_kpath.cpp
)

target_link_libraries(vasp_kpath_tests PRIVATE vasp_kpath GTest::gtest_main)

gtest_discover_tests(vasp_kpath_tests)

add_executable(vasp_surface_tests
    tests/test_surface.cpp
)
target_link_libraries(vasp_surface_tests PRIVATE vasp_surface_lib GTest::gtest_main)
target_compile_definitions(vasp_surface_tests PRIVATE
    TEST_DATA_DIR="${PROJECT_SOURCE_DIR}"
)
gtest_discover_tests(vasp_surface_tests)

add_executable(vasp_eigenval_tests
    tests/test_vasp_eigenval2bands.cpp
    tests/test_vasp_eigenval2bands_i_o.cpp
)

target_link_libraries(vasp_eigenval_tests PRIVATE eigenval2bands GTest::gtest_main)

target_compile_definitions(vasp_eigenval_tests PRIVATE
    TEST_DATA_DIR="${PROJECT_SOURCE_DIR}"
)

gtest_discover_tests(vasp_eigenval_tests)

add_executable(vasp_elastic_tests
    tests/test_elastic.cpp
)
target_link_libraries(vasp_elastic_tests PRIVATE vasp_elastic_lib GTest::gtest_main)
gtest_discover_tests(vasp_elastic_tests)

add_executable(vasp_symmetry_tests
    tests/test_symmetry.cpp
)
target_link_libraries(vasp_symmetry_tests PRIVATE vasp_spglib GTest::gtest_main)
target_compile_definitions(vasp_symmetry_tests PRIVATE
    TEST_DATA_DIR="${PROJECT_SOURCE_DIR}"
)
gtest_discover_tests(vasp_symmetry_tests)

add_executable(vasp_random_tests
    tests/test_random_utility.cpp
)
target_link_libraries(vasp_random_tests PRIVATE vasp_core GTest::gtest_main)
gtest_discover_tests(vasp_random_tests)

# ===== Documentation (Doxygen) =====
find_package(Doxygen OPTIONAL_COMPONENTS dot)
if(DOXYGEN_FOUND)
    set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                   ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generating HTML documentation with Doxygen"
        VERBATIM
    )
    message(STATUS "Doxygen found — build the 'docs' target to generate documentation")
else()
    message(STATUS "Doxygen not found — 'docs' target unavailable")
endif()
