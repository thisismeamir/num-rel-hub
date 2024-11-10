# SpECTRE Installation Guide For AEI Computer (For Usage on Normal Computers)
*A simplified version of the "Setting Up SpECTRE" report for installation.*

## Process:

### Step 1: Install Dependencies

First, make sure that we have the dependencies necessary for the installation process:

```bash
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
  liblapack-dev
```

### Step 2: Pull the Docker Image and Export Files

Next, pull the Docker image and export the files. This method avoids the instability and errors encountered with the GitHub repository versions.

```bash
# Pull the Docker image
docker pull sxscollaboration/spectre:latest

# Create a temporary container
docker create --name temp_container sxscollaboration/spectre:latest

# Copy the SpECTRE directory
docker cp temp_container:/work/ ./your/chosen/path

# Clean up the container
docker rm temp_container
```

Step 2 will provide three folders: `charm_7_0_0`, `spectre`, and `texlive`.

### Step 3: Build Charm++

Navigate into the `charm_7_0_0` folder and run the following commands to build Charm++ for multicore systems:

```bash
# Build Charm++ (for multicore systems)
./build charm++ multicore-linux-x86_64 --with-production --build-shared 

# Set Charm root path (optional)
export CHARM_ROOT=$PWD/multicore-linux-x86_64
```

### Step 4: Build SpECTRE

Next, go into the `spectre` folder and set up the build:

```bash
cd spectre
rm -rf build
mkdir build
cd build

# Configure with CMake
cmake \
  -D CMAKE_C_COMPILER=gcc \
  -D CMAKE_CXX_COMPILER=clang++ \
  -D CMAKE_Fortran_COMPILER=gfortran \
  -D CMAKE_BUILD_TYPE=Release \
  -D CHARM_ROOT=$CHARM_ROOT \
  -D SPECTRE_FETCH_MISSING_DEPS=ON \
  -D MEMORY_ALLOCATOR=SYSTEM \
  ..
```

If this runs without any errors, you're set. Be cautious of warnings, as they can lead to issues with the executables.

### Step 5: Build Components

Now, build the necessary components:

```bash
# Build the CLI interface
make -j$(nproc) cli

# Build all Python bindings
make -j$(nproc) all-pybindings

# Build unit tests
make -j$(nproc) unit-tests

# Build all test executables
make -j$(nproc) test-executables
```

## Important Notes

- **Memory Issues**: Reduce the number of parallel jobs (`-j`) if you experience memory issues during compilation. Even with `-j1`, some executables may encounter memory leaks.
- **Paths**: Adjust paths based on your system configuration.
- **HPC Systems**: For HPC systems, consider using the MPI version of Charm++ (e.g., `mpi-linux-x86_64`) instead of multicore.
- **Debug Mode**: Set `CMAKE_BUILD_TYPE` to "Debug" instead of "Release" for development.
- **Additional CMake Options**: Add options as needed, such as `ENABLE_OPENMP=ON` for OpenMP support.
```