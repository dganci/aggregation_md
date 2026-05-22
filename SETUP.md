## Environments

```bash
conda create -n cg_md python=3.10
conda activate cg_md

cd aggregation_md/envs

conda env update -n unified -f plu_meta_env.yml
conda env update -n unified -f cg_env.yml
conda env update -n unified -f mlcolvar.yml
```

## PLUMED w/ Metatomic

```bash
cd ~

conda activate cg_md

git clone https://github.com/plumed/plumed2.git
cd plumed2

conda install -c conda-forge libmetatomic-torch
```

Add to `~/.bashrc`:
```bash
export LIBTORCH=$HOME/libtorch
export CPATH=$LIBTORCH/include:$LIBTORCH/include/torch/csrc/api/include:$CPATH
export INCLUDE=$LIBTORCH/include:$LIBTORCH/include/torch/csrc/api/include:$INCLUDE
export LIBRARY_PATH=$LIBTORCH/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=$LIBTORCH/lib:$LD_LIBRARY_PATH
if [ -n "$CONDA_PREFIX" ]; then
	export CPATH=$CONDA_PREFIX/include:$CONDA_PREFIX/include/torch/csrc/api/include:$CPATH
 	export LIBRARY_PATH=$CONDA_PREFIX/lib:$LIBRARY_PATH
 	export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
fi
```
Activate: `source ~/.bashrc`

```bash
CPPFLAGS="-I$CONDA_PREFIX/include $CPPFLAGS"
CPPFLAGS="-I$CONDA_PREFIX/include/torch/csrc/api/include $CPPFLAGS"
 
LDFLAGS="-L$CONDA_PREFIX/lib $LDFLAGS"
LDFLAGS="-Wl,-rpath,$CONDA_PREFIX/lib $LDFLAGS"

./configure CC=mpicc CXX=mpicxx CFLAGS="-O2 -fPIC" CXXFLAGS="-O2 -std=c++17 -fPIC -Wno-error -Wno-error=uninitialized -Wno-uninitialized -Wno-unknown-pragmas" MPICC=mpicc MPICXX=mpicxx --enable-mpi --enable-libtorch --enable-libmetatomic --enable-modules=+pytorch:+metatomic LDFLAGS="$LDFLAGS" CPPFLAGS="$CPPFLAGS" --prefix=$HOME/plumed2/

make -j$(nproc)
make install
```

Add to `~/.bashrc`:
```bash
export PLUMED_KERNEL=$HOME/plumed2/lib/libPlumedKernel.so
export PATH=$HOME/plumed2/bin:$PATH
export LD_LIBRARY_PATH=$HOME/plumed2/lib:$LD_LIBRARY_PATH
```
Activate: `source ~/.bashrc`

## GROMACS w/ PLUMED

```bash
cd ~
wget https://ftp.gromacs.org/gromacs/gromacs-2025.0.tar.gz
tar -xzf gromacs-2025.0.tar.gz
cd gromacs-2025.0
mkdir build
cd build
export PLUMED_KERNEL=$HOME/plumed2/lib/libPlumedKernel.so
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/gromacs-plumed -DGMX_MPI=ON -DGMX_OPENMP=ON -DGMX_GPU=OFF -DGMX_BUILD_OWN_FFTW=ON -DGMX_DOUBLE=OFF -DGMX_SIMD=AUTO -DGMX_USE_PLUMED=ON -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DRPATH_POLICY=INSTALL -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
make install
```

Add to `~/.bashrc`:
```bash
source $HOME/gromacs-plumed/bin/GMXRC
```
Activate: `source ~/.bashrc`

## Verify

```bash
which plumed
plumed info --version

which gmx_mpi
gmx_mpi mdrun -h | grep -i plumed
gmx_mpi -version
```
