# Image reconstruction by multivariate kernel density estimation
The application, written in [R language](https://cran.r-project.org/), is devoted to image reconstruction for J-PET,
using [multivariate kernel density estimation (KDE)](https://en.wikipedia.org/wiki/Multivariate_kernel_density_estimation "Wikipedia").

## Prerequisites
Tested on Ubuntu 16.04 LTE with R version 3.4.4 installed. The script operates with a number of shell commands 
(for example, ```readlink -f <incomplete_path_to_file>```), but. basically, all of them are standard.

Additional R packages are required:
* [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html)
* [data.table](https://github.com/Rdatatable/data.table/wiki)
* [ks](https://cran.r-project.org/web/packages/ks/)
* [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html) (for ASCII format only)
* [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html) (optional - for plotting images)

## Input data and parameters 
The application proceeds with [GOJA output format](https://github.com/JPETTomography/j-pet-gate-tools/tree/master/goja#goja-output) 
(16-column ASCII) as the input. The parameters of the reconstruction, as well as the geometry of the scanner, 
are stored in JSON format: it is self-describing and easy to undersand. Please refer to the ```example/``` directory 
for the details.

## Script architecture
The script is designed to work as quick as possible, hence, for data processing, R environments are utilised instead of 
classes:
* ```SetMappingEnvironment.R```
* ```SetAnhPtsEnvironment.R```
* ```SetKDEEnvironment.R```

In comparison, the performance is not crucial for the read of input parameters, so Reference classes are used for 
controllers (wrapped into a list by ```InitParamsFromJSON.R```):
* ```JSONReader.R```
* ```KDEController.R```

More information on the performance of R classes and environments could be found 
[here](https://cran.r-project.org/web/packages/R6/vignettes/Performance.html).

## Usage
There are three main executables: 
* ```MapHitsToCentresOfStrips.R``` -  remaps all XY-coordinates of hits to the centres of strips and saves the result 
to file, according to the format chosen (ASCII or [.RData](http://rfaqs.com/r-workspace-object-image-file)). Only first 
8 columns are saved (pairs of coordinates and times of hits)
* ```DetectAnnihilationPoints.R``` - remaps XY-coordinates of hits to the centres of strips and estimates exact 
positions of annihilation points, using time-of-flight (estimated from times of hits). The result - three columns for 
(X,Y,Z) - is stored according to the format chosen (ASCII or [.RData](http://rfaqs.com/r-workspace-object-image-file)).
* ```EngageKDEReconstruction.R``` - executes KDE on the data for annihilation points. Depends strongly on input JSON file 
and requires ```DetectAnnihilationPoints.R``` to run first.

A single shared argument, a ```.json``` file with parameters is passed to all executables (see the ```example/``` directory):
```
$ Rscript [--vanilla] MapHitsToCentresOfStrips.R <parameters.json>
$ Rscript [--vanilla] DetectAnnihilationPoints.R <parameters.json>
$ Rscript [--vanilla] EngageKDEReconstruction.R <parameters.json>
```
The option ```--vanilla``` [prevents Rscript from
reading R history, profile, or environment files, as well as reloading data or objects from previous sessions](https://stat.ethz.ch/R-manual/R-devel/library/base/html/Startup.html).
