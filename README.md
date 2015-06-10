Follow the installation instructions below to load the HIVBackCalc R package. The package includes [a vignette](https://github.com/hivbackcalc/package1.0/wiki/Vignette-guide-to-the-R-package) that provides a guided analysis. 

## Installing the package

The package is hosted on GitHub and can be installed using the [devtools package](https://github.com/hadley/devtools) by typing the following in R:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github('hivbackcalc/package1.0/HIVBackCalc', build_vignettes=TRUE)
```

## Loading the package

To load the package, type the following in R:

``` r
library(HIVBackCalc)
```

## Vignette

The vignette guides you through a sample analysis. You may view it using the [vignette tab](https://github.com/hivbackcalc/package1.0/wiki/Vignette-guide-to-the-R-package) in this Wiki or from the loaded package by typing

``` r
vignette('Introduction', package='HIVBackCalc')
```
