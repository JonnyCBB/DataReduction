### Next steps:

-	Read sequence file information.
-	Update sigma for observations whose calculated fractions are not (close to) 1
-	This can only be done once sequence information has been read
-	write tests for input parser/Create structure for Reflection/SpaceGroup/Resbin types.

### Inputs still to add:

1) Scale factor from Aimless  
2) Sequence file  
3) Initial structure factor amplitudes from CTruncate

### Known Bugs

-	Running the program with a subset of reflections from a MOSFLM MTZ file that covers more than 1 image will throw an error because it takes the image array that is created doesn't take into account the images from the reflections that are being included. Rather it uses information from ALL batch numbers. This needs to be fixed if someone wants to only use a subset of reflections from their MTZ file. The current alternative would be to run the MTZ file through POINTLESS first and then use the MTZ file generated from the POINTLESS run.

### Things that have to be done at some point.

-	Sort out expected intensities when sequence files aren't given.
-	Sort out the atomic composition when NCS is present.  

### Things that could be done at some point

-	Improve partiality fraction estimation when dealing with partially observed reflections.
-	Could try to come up with a clever way to scale up partially observed intensities that have negatively observed values.
-	Sort out how reflection columns are read from MTZ Dump.
