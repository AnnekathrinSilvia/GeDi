# GeDi 0.99.5

* This version reflects further changes performed upon the Bioconductor reviewing process.
* Removed icons from the vignette to reduce the overall size of the package. 

# GeDi 0.99.4

* This version reflects the changes performed upon the Bioconductor reviewing process
* Added `col_name_genesets` and `col_name_genes` as parameter to the `GeDi()` main app to allow users to specify the relevant column names upon executing the command
* Changes in the R code to comply to best practices (replacing single `|` with `||` and similar)
* Loading the example file does not require anymore the setting of `globalVariables()`
* All files retrieved do use some form of caching for avoiding unnecessary re-download operations
* Reworked the allocation of vectors before `for` loops to avoid unhealthy growing of vectors/matrices

# GeDi 0.99.1

* The handling of the parallelization for the distance calculations is now unified
under the umbrella of BiocParallal, and defaults now to using `SerialParam()` to
avoid unexpected behaviors

# GeDi 0.99.0

* Ready for the submission to Bioconductor!

# GeDi 0.90.0

* Final touches and bug fixes to the main functionality
* Deployment of the package website via pkgdown

# GeDi 0.1.0

* Officially entering the path of the GeDi!
