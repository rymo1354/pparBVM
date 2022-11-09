# pparBVM
Parallelized Bond Valence Method (BVM) parameterization

## Overview
BVM is a localized chemical bonding model that describes the coordination chemistries of inorganic compounds. The BVM Global Instability Index (GII) is calculated from ideal bond length parameters, R0, and bond softness parameters, B, unique to cation-anion pairs. 

This package derives R0 and B from DFT-optimized structures and energies using the framework first presented [here](https://pubs.acs.org/doi/10.1021/acs.jctc.1c01113), whereby the GIIs of DFT ground state structures are minimized while imposing a linear constraint on the relationship between GIIs and DFT energies.

## Dependencies
- Pymatgen
- Numpy
- Scipy
- [pyOpt](https://github.com/madebr/pyOpt)

Additional dependencies listed on the pyOpt Github page.

## pyOpt Installation

After cloning from the Github page, I compiled pyOpt on the [Eagle](https://www.nrel.gov/hpc/eagle-system.html) supercomputing system using the syntax:

``` CC=gcc python setup.py config --f77exec=<PATH> --f90exec=<PATH> build_ext --inplace ```

where `<PATH>` is the path(s) to an appropriate Fortran compiler(s). To use this package, you should add pyOpt to the `$PATH` and `$PYTHONPATH` variables in your `.bashrc`. I had compilation errors when using the default `CC=icc`, so I forced `CC=gcc` to use the GNU Compiler Collection.

## Package Description

### `run_parameterization.py`
Reads in a json-like file type with the structure `{cmpd: {'structures': [pymatgen.core.structure.Structure.as_dict()], 'energies': [float]}}` as a dictionary and prepares this file for parameterization by: 

1. Converting `pymatgen.core.structure.Structure.as_dict()` into `pymatgen.core.structure.Structure`
2. Finding the nearest neighbors of each structure and assigning them to `site_properties`; default method is pymatgen's `CrystalNN` 
3. Creating the structure and energy dictionary, with the same structure as the json
4. Creating the starting parameter dictionary to be optimized with the structure `{'Cation': [pymatgen.core.periodic_table.specie], 'Anion': [pymatgen.core.periodic_table.specie], 'R0': [float], 'B': [float]}`

The structure and energy dictionary and starting parameter dictionaries are then passed to the `BVMParameterizer()` class of `pparBVM/parameterization.py` for parameterization. During parameterization, GIIs are computed using the `GIICalculator()` class of `pparBVM/calculator.py`.

### `submit.py`
Used to submit `run_parameterization.py` to the Eagle computing cluster, which uses [Slurm](https://slurm.schedmd.com/quickstart.html) for job scheduling and management. Can specify the allocation, nodes, cores, etc. as command line arguments. 

## To-Do

1. Implement pyOpt warm restart capabilities
2. Functionalize GIICalculator (anion-anion interactions, second coordination sphere, etc.)
3. Scalability and pyOpt supported optimizer testing 

