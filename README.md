# MORPHEUS
###### Multiparametric Morphometric Analysis of EUcaryotic cellS

## MORPHEUS' quick start guide
###### morphversion 1.0/2019

MORPHEUS plugin can be easily integrated into Fiji menu, simply by copying
`Morpheus_.ijm` file into the Fiji's `plugins` subfolder.

Once installed, the user can launch MORPHEUS from the `Plugins` menu and
MORPHEUS main window should appear.

Through this dialog window, it is possible to specify the following parameters:
* the path of the input directory (i.e., the folder containing the images to be
analyzed);
* the path of the output directory (i.e., the folder in which MORPHEUS will save
the results of the analysis);
* image file suffix (i.e., the extension of the images to be analyzed);
* nucleus file identifier (i.e., a string of characters that must be present in
the file name of all the images stained for nucleus, but NOT in the images
stained for cytoskeleton, allowing MORPHEUS to distinguish among them);
* the anti-spot lower bound epsilon (i.e., the minimum number of pixels a
segmented object must have to not be discarded as 'spot noise');
* tolerance value *T* (i.e., the number of standard deviations from the mean of
the sampling distributions used by MORPHEUS to define the characteristic area
and circularity of the 'typical' isolated cell. The higher the value of T, the
more tolerant the selection process will be);
* 4 checkboxes that allow the user to choose which MORPHEUS' optional functions
to enable among Orientation analysis, Nucleus analysis, Full output and Verbose
log.

**NOTE:** Orientation analysis requires the `OrientationJ` plugin (by D. Sage,
Biomedical Image Group (BIG), EPFL, Switzerland, version 2.0.2 or above) to be
installed. `OrientationJ` latest version can be freely downloaded at
http://bigwww.epfl.ch/demo/orientation/.

When clicking OK button, MORPHEUS starts scanning the user-defined input folder
searching for the images stained for cytoskeleton and no further intervention by
the user is required until the analysis of the whole dataset is completed.

In any case, after MORPHEUS terminates, it is always recommended to check the
results of both segmentation and cell detection processes.
Badly thresholded images can be easily spotted by a quick visual inspection of
the output folder content and their possible presence usually indicates an
insufficient image quality.

If everything worked properly, the message *All output files have been correctly
saved to \<output\>* is printed in Fiji Log window.

For a more detailed guide to MORPHEUS algorithm and output, please refer to:

Ruffinatti FA, Genova T, Mussano F, Munaron L. **MORPHEUS: An automated tool for
unbiased and reproducible cell morphometry**. *J Cell Physiol.*
2020;235:10110–10115. https://doi.org/10.1002/jcp.29768

and **MORPHEUS - User Guide and SM.pdf'**.
