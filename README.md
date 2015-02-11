# package1.0
v1.0 of the R package

## Installation instructions

Download the .tar.gz (Mac OS X, Linux) or .zip file into a directory. Open R in that directory and type

install.packages('HIVBackCalc_1.02.zip', repos=NULL)

(Use the appropriate zipped file in the line above - insert .tar.gz if not on Windows.)

On at least Mac OS X, you may need to include the argument type='source', e.g.

install.packages('HIVBackCalc_1.02.tar.gz', repos=NULL, type='source')

To view the vignette:

library(HIVBackCalc)
vignette('Introduction', package='HIVBackCalc')
