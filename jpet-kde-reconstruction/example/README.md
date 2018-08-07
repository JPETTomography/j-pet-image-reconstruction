# Reconstruction example
This example contains files, essential for the KDE-reconstruction.  In ```data/``` directory, an examplary simulation is placed - for the 1-mm spherical source of the activity 370 kBq, put at (X, Y, Z) = (10, 0, 18.75) [in centimetres] inside 3-layer scanner which reflects the [geometry of the current J-PET prototype](http://koza.if.uj.edu.pl/petwiki/index.php/Big_barrel_dimensions). Z-coordinates and times of hits were smeared according to SiPM readout.

## Input parameters (JSON)
The parameters passed to the script are set in a separate ```.json``` file (```KDE_parameters.json``` in this example). It is composed of three groups of parameters, described below.

#### INPUT-OUTPUT
* **"input-data-path"** - sets the path to the ASCII data file of  [GOJA output format](https://github.com/JPETTomography/j-pet-gate-tools/tree/master/goja#goja-output).
* **"output-file-name"** - assigns the name for the output image (WARNING: the new file will overwrite itself if the name is not changed).
* **"scanner-geometry-file"** - a path to the JSON-file, defining the geometry of the scanner.

#### HIT-MAPPER-OPTIONS
* **"map-from"** - denotes how the hit coordinates will be remapped to the central axis of strips. Allowed options: "*strips_id*" (faster, uses columns 9-10 of [GOJA](https://github.com/JPETTomography/j-pet-gate-tools/tree/master/goja#goja-output)), "*hits*" (slower, directly from Cartesians), "*none*", default is "*strips_id*").
* **"filter-true-events"** - logical, whether to use all coincidences or true only.
* **"save-as-ASCII"** - logical, sets the output format for remapped data (on executing ```MapHitsToCentresOfStrips.R```) and annihilation points (using ```DetectAnnihilationPoints.R```): ASCII if true, .RData otherwise.

#### KS-PACKAGE-OPTIONS
These options are set for multivariate KDE from 'ks' package.
* **"optimise-H-matrix"** - logical, if true, the optimised value (estimated experimentally) is set for plug-in bandwidth selector (H-matrix), otherwise the default one - see [the topic on Hpi in documentation for 'ks'-package](https://cran.r-project.org/web/packages/ks/ks.pdf).
* **"normalise-by"** - how to normalise the intensity of the image (from 0 to 1). Allowed options: "*peak*" (1 for the brightest pixel), "*sum*" (the sum of all pixel intensities equals to 1).
* **"resolution"** - defines the resolution of the resulting 3D image.
* **"axis-ranges"** - allows KDE being run on sub-region. If "*auto-detect*" is set to false, the ranges for Cartesians are to be defined, otherwise ("*auto-detect: true*") the whole field-of-view is in count.
* **"save-as-ASCII"** - logical, if true, the reconstructed image will be stored in a separate directory, comprising ASCII files for all axis (coordinates of pixels) and 3D image (as a table **[X_pixel_ID, Y_pixel_ID, Z_pixel_ID, Pixel_Intensity]**). If false, the data will be stored as R workspace ```.RData```.

## Scanners (JSON)
The directory ```scanner/``` contains two ```.json``` files which define geometries, for the ideal 384-strip scanner and for the current J-PET prototype (big barrel). The latter is required for the example (```data/LAB_192_10_00_19_SI```).

**(!) Important note!** The script is consistent with the cylindrical geometry only and would need serious redesign in order to match other perspective configurations of J-PET scanner.

The basic parameters of the scanner ```.json``` file are self-describing and easy to understand. The only puzzling feature is ```angle-bias-coefficient```: it defines the angle between the first detector/strip in a layer and X-axis, calculated as *angle-bias-coefficient \* pi / no_of_strips_in_layer*. For the big barrel, these angles are 0, 3.75 and 1.875 degrees (for 1st, 2nd and 3rd layer, respectively). See more on this on PET wiki: [here](http://koza.if.uj.edu.pl/petwiki/images/5/5b/Plyta_katy.CATDrawing.pdf) and [here](http://koza.if.uj.edu.pl/files/46d8a923615abf37abff1c0e1409137e/oPs_jpet.pdf).

## Running the example on CIŚ cluster
The sequence of steps below is shown for the interactive mode (```qsub -I```), but could also be wrapped up into ```.sh``` script to be executed as a queue.
```
[user@usrint2 ~]$ qsub -I
[user@wnXXXXX ~]$ cd jpet-kde-reconstruction/ # if cloned from github
...$ module load R/3.4.1
...$ Rscript --vanilla src/DetectAnnihilationPoints.R example/KDE_parameters.json
...$ Rscript --vanilla src/EngageKDEReconstruction.R example/KDE_parameters.json
```

## Plotting images
The example for plotting cross-cections (XY and XZ) of 3D image could be found in the file ```ExampleForPlottingImages.R```. Unlike other scripts, it shoud be run interactively from R and not from the shell. In case of using the cluster:
```
$ cd jpet-kde-reconstruction/ # if cloned from github
$ module load R/3.4.1
$ R
R version 3.4.1 (2017-06-30) ...
...
> setwd("example/")
> load("ReconstructedImage.RData")
> str(image_object)
> image_object[["H"]] # H-matrix
...
> q() # quits R
```
The plots (if made) will be stored in a file ```Rplots.pdf```. There's also an alternative of plotting images from within R using [pdf, png, dev.off](http://rfunction.com/archives/812).

**(!) Impornant note!** Standard R functions does not support easy plot with the legend (for colour scale), so it's better to use ['ggplot2' package](http://r-statistics.co/ggplot2-Tutorial-With-R.html).
