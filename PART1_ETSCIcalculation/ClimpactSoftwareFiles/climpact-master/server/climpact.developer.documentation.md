# Climpact developer document

This document provides a high level description of the [Climpact](https://github.com/ARCCSS-extremes/climpact) source code as of November 2017. 

Climpact is an R based program that runs on the user's computer, primarily via the web browser. It uses the R [Shiny](http://shiny.rstudio.com) package to do this. Climpact is also dependent on numerous other R packages which are installed by the installer script (see the README for instructions).

While R Shiny is a web-app platform which makes Climpact run through a browser window, Climpact does not require an internet connection to operate once it is installed.

## Versioning

Version numbering follows the tpyical X.Y.Z pathway where X denotes a major change in functionality and/or appearance, Y a minor change in functionality and/or appearance and Z a bug fix. Version number should be updated in climpact.etsci-functions.r whenever changes are pushed to github.

## Source files

All files in the Climpact home directory (and its subdirectories) ending in .r and .R constitute Climpact source files. Due to the design philosophy of R Shiny, Climpact source files are largely divided into two categories; UI (user interface) and server files. The former governs the rendering of the Climpact interface and the latter contains the core Climpact functionality. These files are stored in the ui/ and server/ subdirectories, respectively.

Initiating Climpact, via the runApp() command (again, see the README for instructions), accesses the app.R file in the Climpact home directory, which links to other UI files in the ui/ subdirectory. There are other .r files in the Climpact home directory. Most of these are *.wrapper.r scripts, which allow access to Climpact's more advanced functions (namely gridded indices and batch processing) via the Linux/MacOs/Windows command line. This functionality is also available - though not 100% user proof - via the user interface.

## Making changes to Climpact

Due to the increasing complexity of Climpact, modifications to the software require a reasonable time to gain familiarity. Particularly, for users not familiar with R Shiny and wanting to change aspects of the front-end will require some time to appreciate the programming philosophy that R Shiny follows ('reactive' programming).

Needless to say, modifications to the user interface occur via the app.R file and the Climpact source files in the ui/ subdirectory, and requires a basic understanding of R Shiny. It is recommended to watch the [R Shiny tutorial](https://shiny.rstudio.com/tutorial/).

### Adding new indices

Adding new indices to Climpact should not require knowledge of R Shiny. However, it does require knowledge of R packages climdex.pcic and - particularly - climdex.pcic.ncdf. Climpact makes use of both of these packages for calculating the indices and the latter is reasonably complex.

While Climpact depends on the publically available version of climdex.pcic, it relies on a modified version of climdex.pcic.ncdf. This is because additional indices included in Climpact but not included in the base climdex.pcic package are simply added via the climpact.GUI-functions.r file, which calls necessary functions from climdex.pcic when needed. However, for calculating new indices on netCDF files the climdex.pcic.ncdf package had to be modified. Thus, when Climpact is installed it installs a modified version of climdex.pcic.ncdf and overwrites any previously installed version (with warning). 

How complex the code changes need to be when adding a new index will depend on how complex the new index is to calculate and write out. If it has similar inputs and outputs to current indices then the changes are minimal, however, for example, when the SPEI and heatwave indices were added, significant modifications were needed to pass the needed data to the index functions and then to write out netCDF files with the correct structure.

### Adding a new function for your index

The first step to adding a new index to Climpact is to create the function that calculates the index. This should be placed inside server/climpact.etsci-functions.r. See other examples of ET-SCI indices added in this file. New indices also need to be added to server/climate.indices.csv. 

### Modifying climdex.pcic.ncdf

The second step in adding a new index to Climpact requires modifying the climdex.pcic.ncdf package that comes with Climpact. The modified version of climdex.pcic.ncdf that comes with Climpact can be found at server/pcic_packages/climdex.pcic.ncdf.climpact.tar.gz. To add new indices to this package in order that they are calculated on netCDF files follow these general steps:
1. Extract this file to create a climdex.pcic.ncdf folder.
2. Go to climdex.pcic.ncdf/R. Here you will see the file ncdf.R, this is the package's source code. As well you will see a folder named orig/ with a copy of the base version of ncdf.R released by PCIC.
3. When adding indices to Climpact it is instructive to look at the difference between the above source files (e.g. vimdiff) to see how ET-SCI indices were added to the original code.
4. Code changes to ncdf.R. A summary of the changes required for simple indices are described below (excuse any error or vagueness here due to the author's memory of this process):
* In get.climdex.variable.list() you will see lists that need to be edited to identify which indices can be calculated with which of the three input variables TX, TN and/or PR, as well as which indices are calculated monthly-only and annual-only. Add your index to the appropriate lists.
* In get.climdex.functions() default parameter values for your index's function can be added. You will also see list() variables defined for indices that require parameters to be sent to each 'worker' (in a multicore operation).
* Once you modify the source file you will want to remove the currently installed climdex.pcic.ncdf package, create a new package file and install that. Follow these steps;
..1. Create a new tar ball of your modified code by running *tar -czvf climdex.pcic.ncdf.climpact.tar.gz climdex.pcic.ncdf/*
..2. Open R and run *remove.packages("climdex.pcic.ncdf") ; install.packages("climdex.pcic.ncdf.climpact.tar.gz",repos = NULL, type="source")*
5. For complex indices (e.g. that have multiple output fields or require new thresholds to be calculated) it is recommended that the changes made for the SPEI and heatwave indices are studied in detail (e.g. vimdiff).

## TO-DO list for future updates

* Improve exception handling when errors occur for gridded indices and batch processing. This should include displaying *realtime* console output to the browser while indices are being calculated. Current version of Shiny inhibits this. **A decent step in the right direction is to remove the need for a system call to 'Rscript' when calculating gridded indices, as this prevents proper errors messages being returned to the R console**.
* Allow indices to be 'loaded' into Climpact when calculating correlations with sector data. Currently, the *Load and Check Data* and *Calculate Indices* processes have to be completed before the user can calculate correlations with sector data, even if the indices have been calculated previously. An alternate and possibly better option is to have the ability to load sector data in the Load and Check Data tab, and thus have correlations calculated automatically from the beginning of the process.
* Allow monthly sector data to be loaded and correlated with indices.
* Implement simple regression plots between sector data and indices.
* New index: % of mean precipitaiton. Similar to Tosi's index.
* For future indices requiring humidity, utilise published proxies of relative humidity based on temperature and precipitation.
* **clean code: reduce usage of global variables, improve functionalisation of the code.**
