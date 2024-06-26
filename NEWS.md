# simdata v0.4.0
Release for CRAN with updates in vignettes to reflect recent changes to nhanesA
package.

# simdata v0.3.0.9004
Adding automated quantile function estimation for NORTA via the 
`quantile_functions_from_data()` utility.

Fixed errors about missing `colapply_functions()` in specific situations.

# simdata v0.3.0.9003
Implemented partial functions.

# simdata v0.3.0.9002
Updating defaults and functionality for correlation network visualisation.

# simdata v0.3.0.9001
Added tolerance to correlation matrix check.

# simdata v0.3.0.9000 (build date 2021-04-09)
Added NORTA based simulation design and associated functionality. Added
vignettes explaining the workflow, and the technical design of the package.
 
Additions: 

- simdesign_norta to simulate data based on a NORTA design
- optimize_cor_mat to facilitate NORTA simulation
- Vignette demonstrating the NORTA workflow
- Vignette explaining the technical implementation of the package

Changes: 

- Updated vignettes
- Changes in function names to improve naming consistency and enhance 
    searchability: 
    - mvtnorm_simdesign -> simdesign_mvtnorm 
    - discunif_simdesign -> simdesign_discunif
    - names_from_function_list -> get_names_from_function_list

# simdata v0.2.0.9003 (build date 2021-04-08)
Changes:

- Post-processing pipeline: renamed process_data function to 
    do_process; renamed process_truncate to process_truncate_by_iqr
    
Additions:

- discunif_simdesign to simulate circular data
- process_truncate_by_threshold function to truncate by fixed thresholds
- names_from_function_list to obtain names of functions from function_list
- is_cor_matrix to check if a matrix is a correlation matrix

Fixes and enhancements:

- More robust error handling in simdesign
- More robust and transparent naming of variables
- Basic error checking in essential functions
- More robust handling of function_lists

# simdata v0.2.0.9002 (build date 2020-05-14)
Renamed the package to simdata. 

Fixes:

- Fixed cor_from_upper to work with single vectors

# simdata v0.2.0 (build date 2019-11-26)
Numerous fixes and enhancements. If you have worked with previous versions of 
the package it is likely that you will have to adapt your code to the new
interfaces implemented in this version. However, the necessary changes should be 
minimal.

Changes:

- Reworked simdesign class
- Reworked simulate_data and simulate_data_conditional function
- Reworked network visualization

Additions: 

- Package Demo vignette to demonstrate functionality

# simdata v0.1.0 (build date 2019-11-14)
Initial release (named simulatoR).