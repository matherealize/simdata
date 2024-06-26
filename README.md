# `simdata`: An R package to create simulated datasets

This package provides lightweight and easy-to-use data generating mechanisms
for statistical simulation studies.

## Installation

``` r
install.packages("simdata")

# Alternatively, install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("matherealize/simdata")

# to also build vignettes locally use
# devtools::install_github("matherealize/simdata", build_vignettes = TRUE)
```

In case you encounter an http 401 error during installation from GitHub, this is
not due to the package itself, but due to rate limits set by GitHub. Please 
refer to [these steps from the remotes Github repository](https://github.com/r-lib/remotes/issues/330#issuecomment-578474009) 
to resolve the issue.

## Usage

Please take a look at the package vignettes:

- [Introduction to `simdata` package and basic workflow](https://matherealize.github.io/simdata/articles/Demo.html)
- [Introduction to NORTA based simulation workflow](https://matherealize.github.io/simdata/articles/NORTA_demo.html)
- [Technical documentation, implementation overview](https://matherealize.github.io/simdata/articles/Technical_documentation.html)

## Status and future plans

The package is in active development, but is already used in several simulation 
studies. We plan to provide a library of data generators for interesting 
real-world dataset in upcoming versions of the package.

## Contact

For feedback and issues when working with the package please create a new
[Issue at Github](https://github.com/matherealize/simdata/issues) or 
send [me an email](mailto:michael.kammer@meduniwien.ac.at).
