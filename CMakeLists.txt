cmake_minimum_required (VERSION 3.2)

project (hw5 LANGUAGES CXX) 

# Compiler requirements
set(CMAKE_CXX_STANDARD 20) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

# Find JSON Libraries
#find_package(nlohmann_json REQUIRED)

# FIND HDF5 Libraries 
#find_package(HDF5 REQUIRED)
#find_package(HighFive REQUIRED)

# Find Linear Algebra Header Libraries
find_package(Armadillo REQUIRED)
# find_package(Eigen3 3.4 REQUIRED NO_MODULE)

# Find Linear Algebra Runtime Libraries
# find_package(BLAS REQUIRED)
# find_package(LAPACK REQUIRED)

# Includes, Links and Defines 
# link_libraries(${ARMADILLO_LIBRARIES})
# link_libraries(Eigen3::Eigen) 
# link_libraries(armadillo)
# link_directories(nlohmann_json::nlohmann_json)
# link_libraries(HDF5::HDF5)
# add_compile_definitions(ARMA_USE_HDF5)

# add subdir 
#add_subdirectory(src)

# Add your existing include directories
include_directories(${PROJECT_SOURCE_DIR}/Include)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Add libraries built from src/*.cpp
add_library(Atom ./src/Atom.cpp)
add_library(AtomReader ./src/AtomReader.cpp)
add_library(Molecule_init ./src/Molecule_init.cpp)
add_library(Molecule_utils ./src/Molecule_utils.cpp)
add_library(Molecule_overlap ./src/Molecule_overlap.cpp)
add_library(Molecule_HF ./src/Molecule_HF.cpp)
add_library(Molecule_AnalytGrad_Ecndo ./src/Molecule_AnalytGrad_Ecndo.cpp)


# Link librarues, Ensure Armadillo is linked with these targets
target_link_libraries(Molecule_AnalytGrad_Ecndo PRIVATE 
                                ${ARMADILLO_LIBRARIES}
                                Molecule_init 
                                Molecule_utils
                                Molecule_overlap
                                Molecule_HF)
target_link_libraries(AtomReader PRIVATE Atom )

# Add executables for tests
add_executable(hw5_2 ./src/hw5_2.cpp)

# Link libraries to tests and include necessary directories
target_include_directories(Atom PRIVATE ./Include)
target_include_directories(AtomReader PRIVATE ./Include)
target_include_directories(Molecule_AnalytGrad_Ecndo PRIVATE ./Include)

# Link executables to their respective libraries and Armadillo
target_link_libraries(hw5_2 PRIVATE AtomReader Molecule_AnalytGrad_Ecndo)