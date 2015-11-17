### Next steps:

-	Write Unscented Kalman Filter algorithm
-	Create Parameters for the Unscented Kalman Filter algorithm
-	write tests for input parser/Create structure for Reflection/SpaceGroup/Resbin types.

### Inputs still to add:

1) Variance coefficients so the user can manually adjust the size of the variances

### Known Bugs

-	Running the program with a subset of reflections from a MOSFLM MTZ file that covers more than 1 image will throw an error because the image array that is created doesn't take into account the images from the reflections that are being included. Rather it uses information from ALL batch numbers. This needs to be fixed if someone wants to only use a subset of reflections from their MTZ file. The current alternative would be to run the MTZ file through POINTLESS first and then use the MTZ file generated from the POINTLESS run.
-	Some reflections can be allocated a centroid phi value that isn't on an image on which it was observed (I don't know how this works but MOSFLM allows it). I have to figure out how to deal with this.

### Things that have to be done at some point.

- I need to sort out the implementation of the Rician distribution for the structure factor amplitudes. The bessel function throws an error (Base.Math.AmosExecption(2) at special/bessel.jl:142) with some of the numbers that are used. I don't know what it is but it seems to be a problem when the Gaussian mean is much larger than the standard deviation. Maybe check it against another programming language to see what the problem might be.
-	Check to make sure the reflection multiplicity (epsilon factor) is correct
-	Sort out expected intensities when sequence files aren't given.
-	Sort out the atomic composition when NCS is present.
-	Outlier rejection
-	Sort out error messages so they use '@sprintf' macros. At the moment they wont give the correct strings in general
-   Sort out image rejection
-   Properly comment the code.
-   Actually output log information to the console so that the user knows what's going on.
- Sort out rounding of values when sorting out the image numbers during the MTZ file parsing. I've left a "*****" comment where this needs to be sorted. 

### Things that could be done at some point

-   Add some optimisation for the variances of the process and the observation matrix.
-   I need to take into account the uncertainties of our estimates of the scale factor and the B factor. Looks like this could be done with application of the total law of variance or marginalising the conditional probabilities.
-	My object structure isn't great. Both Reflection and Diffraction image types contain reflection Observations. This is a waste of memory. I may have to restructure this so that only Diffraction Images contain observations of reflections because that's where it seems sensible.
-	Improve partiality fraction estimation when dealing with partially observed reflections.
-	Could try to come up with a clever way to scale up partially observed intensities that have negatively observed values.
-	Sort out how reflection columns are read from MTZ Dump.
-	Be explicit with reading in the column information from the MTZ Dump output from the CTruncate file. At the moment I've just looked at the column numbers and inserted the numbers straight from the output file. This may not be consistent if other input files are read.
-   Perform a weighted average when calculating the mean intensity of resolution bins from the images.
