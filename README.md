![FMI.jl Logo](https://github.com/JuliaFMI/FMI.jl/blob/main/logo/dark/fmijl_logo_640_320.png "FMI.jl Logo")
# FMIZoo.jl

## What is FMIZoo.jl?
[*FMIZoo.jl*](https://github.com/JuliaFMI/FMIZoo.jl) is a collection of testing and example FMUs ([fmi-standard.org](http://fmi-standard.org/)) for the Julia libraries [*FMI.jl*](https://github.com/JuliaFMI/FMI.jl) and [*FMIFlux.jl*](https://github.com/JuliaFMI/FMIFlux.jl). 

[![Test (latest)](https://github.com/JuliaFMI/FMIZoo.jl/actions/workflows/TestLatest.yml/badge.svg)](https://github.com/JuliaFMI/FMIZoo.jl/actions/workflows/TestLatest.yml) 
[![Test (LTS)](https://github.com/JuliaFMI/FMIZoo.jl/actions/workflows/TestLTS.yml/badge.svg)](https://github.com/JuliaFMI/FMIZoo.jl/actions/workflows/TestLTS.yml)
[![Coverage](https://codecov.io/gh/JuliaFMI/FMIZoo.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaFMI/FMIZoo.jl)

## How can I use FMIZoo.jl?
1. open a Julia-Command-Window, activate your preferred environment
1. goto package manager using ```]```
1. type ```add FMIZoo``` or ```add "https://github.com/JuliaFMI/FMIZoo.jl"```

## What FMI.jl-Library to use?
![FMI.jl Logo](https://github.com/JuliaFMI/FMI.jl/blob/main/docs/src/assets/FMI_JL_family.png "FMI.jl Family")
To keep dependencies nice and clean, the original package [*FMI.jl*](https://github.com/JuliaFMI/FMI.jl) had been split into new packages:
- [*FMI.jl*](https://github.com/JuliaFMI/FMI.jl): High level loading, manipulating, saving or building entire FMUs from scratch
- [*FMIImport.jl*](https://github.com/JuliaFMI/FMIImport.jl): Importing FMUs into Julia
- [*FMIExport.jl*](https://github.com/JuliaFMI/FMIExport.jl): Exporting stand-alone FMUs from Julia Code
- [*FMIBase.jl*](https://github.com/JuliaFMI/FMIBase.jl): Common concepts for import and export of FMUs
- [*FMICore.jl*](https://github.com/JuliaFMI/FMICore.jl): C-code wrapper for the FMI-standard
- [*FMISensitivity.jl*](https://github.com/JuliaFMI/FMISensitivity.jl): Static and dynamic sensitivities over FMUs
- [*FMIBuild.jl*](https://github.com/JuliaFMI/FMIBuild.jl): Compiler/Compilation dependencies for FMIExport.jl
- [*FMIFlux.jl*](https://github.com/JuliaFMI/FMIFlux.jl): Machine Learning with FMUs
- [*FMIZoo.jl*](https://github.com/JuliaFMI/FMIZoo.jl): A collection of testing and example FMUs

## What Platforms are supported?
[*FMIZoo.jl*](https://github.com/JuliaFMI/FMIZoo.jl) is tested (and testing) under Julia Versions *1.6.7 LTS* (64-bit) and *latest* (64-bit) on Windows *latest* (64-bit) and Ubuntu *latest* (64-bit). MacOS and Julia (32-bit) should work, but untested.