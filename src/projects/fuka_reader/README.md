## FUKA Reader

Read FUKA ID for BNS, BBH, and BHNS binaries.
fuka_reader interfaces with [fukaccia](https://github.com/unkaktus/fukaccia),
which performs interpolation using FUKA functions.

### Installation

The easiest way to install fuka_reader without compiling any C++ code
is to install the binary fukaccia package using Mamba.
Alternatively, you can build fukaccia yourself by following the instructions
on its GitHub page, finding all the dependencies, and readjusting library paths in
the instructions below.

To install MambaForge, run:
```shell
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

Most likely, the cluster you are working on does not have internet access,
so you would need to install `mitten` locally; see https://github.com/unkaktus/mitten.
Then, you can enable internet access by mitten-ning into the cluster instead of ssh-ing.
Also, you can then directly clone BAM or download anything from the internet via mitten.

Install fukaccia with all its dependencies:

```shell
mitten cluster
mamba install -c https://mamba.unkaktus.art fukaccia
```

Clone BAM using your access token:
```shell
git clone https://username:token@git.tpi.uni-jena.de/bamdev/BAM.git
```

Add `fuka_reader` project to your `MyConfig` file:
```Makefile
projects += src/projects/fuka_reader
```
and add the `fukaccia` dependency:
``` Makefile
SPECIALLIBS += -lfukaccia -lfuka_exporter -lfftw3 -lgsl -lblas
```

Use your favorite BAM project cloner, for example, the Makefile.

Build BAM:
```shell
make new
```

If you compiling with non-mamba GCC, you need to add locations to your MyConfig:
```
libsys += -lrt
INCS += -I${CONDA_PREFIX}/include
SPECIALLIBS += -L${CONDA_PREFIX}/lib
```
and export the location for the required libraries
before running BAM:

```shell
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib
```

### Usage

There are the following parfile paramters:

* `fuka_reader_binary_type` - Type of the binary [BBH, BNS, BHNS]
* `fuka_reader_info_filename` - Path _to_ the FUKA .info file for the computed ID

Parameters relevant only for BHs:
* `fuka_reader_interpolation_offset` - Interpolation offset (in units of r_AH)
* `fuka_reader_interpolation_order` - Interpolation order for smooth junk
* `fuka_reader_relative_dr_spacing` - Relative dr spacing for the interpolation polynomial



### Author
Ivan Markin, 09/2023
