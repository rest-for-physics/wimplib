[![DOI](https://zenodo.org/badge/324291710.svg)](http://doi.org/10.5281/zenodo.4528985)
[![website](https://img.shields.io/badge/user-guide-E8B6FF.svg)](https://rest-for-physics.github.io)
[![api](https://img.shields.io/badge/user-API-FFCA78.svg)](https://sultan.unizar.es/rest/)
[![forum](https://img.shields.io/badge/user-forum-AAFF90.svg)](https://rest-forum.unizar.es/)

### WIMPlib

This is a REST-for-Physics library for the calculation of different WIMP parameters and sensitivity plots using REST metadata structures.

### RestWIMPLib installation

Once you have all the prerequisites installed you need to add the library at the REST-for-Physics framework compilation stage. Go to the main framework build directory and add the WIMP library as a compilation option.

```
cd framework/build
cmake ../ -DCMAKE_INSTALL_PREFIX=../install/ -DRESTLIB_WIMP=ON
make -j4 install
```

### Publications

This repository makes use of the following published codes:
- K. Altenmuller et al, REST-for-Physics, a ROOT-based framework for event oriented data analysis and combined Monte Carlo response, [Computer Physics Communications 273, April 2022, 108281](https://doi.org/10.1016/j.cpc.2021.108281).
