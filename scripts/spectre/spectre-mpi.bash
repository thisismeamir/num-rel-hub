#!/bin/bash

# Update package list and install dependencies
sudo apt update
sudo apt install -y \
  build-essential \
  clang \
  gfortran \
  cmake \
  libboost-all-dev \
  libgsl-dev \
  libhdf5-dev \
  python3-dev \
  python3-pip \
  libblas-dev \
  liblapack-dev \
  mpich

# Pull SpECTRE Docker image
docker pull sxscollaboration/spectre:latest

# Create and export the SpECTRE files from Docker
docker create --name temp_container sxscollaboration/spectre:latest
docker cp temp_container:/work/ ./spectre_install
docker rm temp_container

# Navigate to Charm++ folder and build for MPI
cd ./spectre_install/charm_7_0_0
./build charm++ mpi-linux-x86_64 --with-production --build-shared

# Set Charm root path
export CHARM_ROOT=$PWD/mpi-linux-x86_64

# Configure SpECTRE with CMake
cd ../spectre
rm -rf build
mkdir build
cd build
cmake \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_Fortran_COMPILER=mpifort \
  -D CMAKE_BUILD_TYPE=Release \
  -D CHARM_ROOT=$CHARM_ROOT \
  -D SPECTRE_FETCH_MISSING_DEPS=ON \
  -D MEMORY_ALLOCATOR=SYSTEM \
  ..

