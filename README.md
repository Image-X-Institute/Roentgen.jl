# Dose Calculations
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://acrf-image-x-institute.github.io/DoseCalculations.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://acrf-image-x-institute.github.io/DoseCalculations.jl/dev/)

*Lars Mejnertsen*

A Julia package for calculating dose for radiotherapy treatments.

This package is under development, so the implementation may change between versions.

![dose_recon_example](docs/src/assets/dose-reconstruction.png)

## Setup and Usage

To add the package to your environment, open the Julia REPL and type the `]` to open the package manager.
Add the package by entering:

```julia
add https://github.com/ACRF-Image-X-Institute/DoseCalculations.jl
```

Once completed, you can exit the package manager by pressing Backspace.

You can then use the package by entering:

```julia
using DoseCalculations
```

There are examples for dose calculations in the `examples/` folder

## Running the Examples

The `examples/` directory contains its own environment to ensure you have all the packages required to run the code.

1. Navigate to `examples/` in your terminal

2. Modify `path/to/data` strings in the example file to point to the required data.

3. Run the following command, swapping out `{example_file}` with the name of the file:

    ```sh
    julia --project=. {example_file}.jl
    ```

**Note:** Julia is a Just-In-Time (JIT) compiled language. This means that every time you run the scripts from the terminal, it re-compiles the code, causing the code to run slower than expected. If you want to run these scripts multiple times, a better way to run this code would be to use the REPL or Visual Studio Code with the Julia Extension.

## Documentation

The documentation is still under development.
The documentation is generated by [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/) and located in the `docs/` folder.

To build:

1. Navigate to the project directory in the terminal.
2. Run `julia --project=docs make.jl`

This will generate the documentation HTML in `docs/build`.

## Directory Structure

```
Project.toml	Dependency Files
src/			Source Files
examples/		Example Usage Scripts
docs/			Documentation
```
