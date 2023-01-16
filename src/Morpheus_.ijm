// MORPHEUS -- release Sep-2020
// ImageJ 1.x macro language version
// OrientationJ 2.0.2 (or above) required

//---------------------------------------------------------------------------------------------------//
// Morpheus GUI
//---------------------------------------------------------------------------------------------------//

// Universal I/O parameter notation (ImageJ2 required)
#@ File(label="Input directory", style="directory") input
#@ File(label="Output directory", style="directory") output
#@ String(label="File suffix", value=".tif") suffix
#@ String(label="Nucleus file identifier", value="dapi") identifier
#@ Integer(label="Anti-spot lower bound (pixel^2)", value=200) lowBound
#@ String(label="Tolerance (sigma)", choices={"1","2","3","4","5","6"}, style="listBox", value="3") tol
#@ String (label="Orientation analysis", choices={"None", "Double-Scale", "All Filaments"}, style="radioButtonHorizontal") orientStr
#@ Boolean(label="Nucleus analysis", value=false) nucleus
#@ Boolean(label="Full output", value=false) full
#@ Boolean(label="Verbose log", value=true) verbose

//---------------------------------------------------------------------------------------------------//
// Preliminaries
//---------------------------------------------------------------------------------------------------//

// Batch Mode! ...to speed up the macro
setBatchMode(true);

// For Reproducibility (when adding noise)
random("seed", 101);

// Close all open images and windows
run("Close All");
if (isOpen("Results")) {
	selectWindow("Results");
    run("Close");
}

// Check ImageJ version:
// ImageJ2 (or Fiji) is required for #@ parameter notation
ver = getVersion(); // ImageJ -> 1.53e, while ImageJ2 -> 2.1.0/1.53e
verArray = split(ver,"./"); // Get an array of strings: 2,1,0,1,53e
verInt = parseInt(verArray[0]); // Convert the first element to an integer
if(verInt < 2)
	exit("Morpheus requires ImageJ2 (or Fiji).");
// 1.52a (or higher) is required for Table Functions
// version 1.52k (or higher) is required for setOption("ScaleConversions", boolean) functions
requires("1.52k");

// Morpheus Version
morphversion = "1.2/Oct-2020";

// Set Default Options
setDefaultOptions(); // A Morpheus' function - see definition below

// Remap orientStr to boolean
if (orientStr == "None")
	orient = false;
else
	orient = true;

// Print General Information as Log
systeminfo(morphversion); // A Morpheus' function - see definition below

// Force case-insensitiveness
suffix = toLowerCase(suffix);
identifier = toLowerCase(identifier);

// Define equivalence classes for ambiguous extensions
suffix2 = suffix;
if(endsWith(suffix, "tiff") || endsWith(suffix, "tif")) {
	suffix = ".tif";
	suffix2 = ".tiff";
}
if(endsWith(suffix, "jpeg") || endsWith(suffix, "jpg")) {
	suffix = ".jpg";
	suffix2 = ".jpeg";
}

// Safety option
if(lengthOf(identifier) == 0)
	identifier = "dapi";

// Organize Output Directory
File.makeDirectory(output + File.separator + "Cell");
if(nucleus)
	File.makeDirectory(output + File.separator + "Nucleus");

// Remap Tolerance from string to integer
tol = parseInt(tol);

// Global variables
var objCounter = 0; // Count all the objects detected
var shortlist = newArray(); // The subset of "list" containing only those files suitable for morphometric analysis
var degreeAxis = newArray(); // -90:1:89 array
for (i = 0; i < 180; i++)
	degreeAxis = Array.concat(degreeAxis, i-90);

//---------------------------------------------------------------------------------------------------//
// Main
//---------------------------------------------------------------------------------------------------//

// Learn the characteristic features of the single cell from the whole dataset
run("Set Measurements...", "area shape redirect=None decimal=3");
// "pop" is an array containing median areas and extreme circularities of each sample image (sampling distributions)
pop = GetSamplingDistributions(input, suffix, suffix2, identifier, lowBound); // A Morpheus' function - see definition below

// Check for contents
if(pop.length == 0) {
	messageOut = "No suitable " + suffix + " files have been detected in the folder\n" + input;
	print(messageOut);
	exit(messageOut);
}

// Separate areas from circularities
sampleNum = pop.length/3;
popAreas = newArray(sampleNum);
popMinC = newArray(sampleNum);
popMaxC = newArray(sampleNum);
for (i = 0; i < sampleNum; i++) {
	popAreas[i] = pop[3*i];
	popMinC[i] = pop[3*i+1];
	popMaxC[i] = pop[3*i+2];
}

// Log the two sampling distributions
if(verbose) {
	print("\nArea - sampling distribution of the sample median:");
	Array.print(popAreas);
	print("\nCircularity - sampling distribution of the sample minimum:");
	Array.print(popMinC);
	print("\nCircularity - sampling distribution of the sample maximum:");
	Array.print(popMaxC);
}

// Define typical cell area and circularity range using the means of their corresponding sampling distributions 
Array.getStatistics(popAreas, minNoUse, maxNoUse, typicalArea, AreaStd);
typicalRadius = sqrt(typicalArea/PI);
Array.getStatistics(popMinC, minNoUse, maxNoUse, typicalMinC, MinCStd);
Array.getStatistics(popMaxC, minNoUse, maxNoUse, typicalMaxC, MaxCStd);
if(verbose) {
	print("\n");
	print("Characteristic Area of the single cell +/- Standard Error = (", typicalArea, " +/- ", AreaStd/sqrt(sampleNum), ") pixel^2", " --- Relative = ", 100*(AreaStd/sqrt(sampleNum))/typicalArea, "%");
	print("Characteristic Radius of the single cell +/- Standard Error = (", typicalRadius, " +/- ", (AreaStd/sqrt(sampleNum))/(2*PI*typicalRadius), ") pixel");
	print("\n");
}

// Define new bounds for single cell size and circularity
/* e.g. if tolerance = 3
 * (typicalArea + 3*Standard_Error) = 99.7% Confidence Interval upper bound for the "mean median" statistic
 * 2*(typicalArea + 3*Standard_Error) is a conservative estimate for the predicted typical size of two-cell complex
 * in other words, if an object has a size > highBound can be legitimately suspected to be a complex of two cells
 */
highBound = 2*(typicalArea + tol*AreaStd/sqrt(sampleNum));
if(0.5*(typicalArea-tol*AreaStd/sqrt(sampleNum)) > lowBound)
	lowBound = 0.5*(typicalArea - tol*AreaStd/sqrt(sampleNum));

// High tolerance (6):	typicalMinC + 0*MinCStd AND typicalMaxC - 0*MaxCStd
// Mid tolerance (3):	typicalMinC + 3*MinCStd AND typicalMaxC - 3*MaxCStd
// Low tolerance (1):	typicalMinC + 5*MinCStd AND typicalMaxC - 5*MaxCStd
lowBoundC = typicalMinC + (6-tol)*MinCStd/sqrt(sampleNum);
highBoundC = typicalMaxC - (6-tol)*MaxCStd/sqrt(sampleNum);
boundFlag = false;
if(highBoundC <= lowBoundC) {
	lowBoundC = typicalMinC;
	highBoundC = typicalMaxC;
	boundFlag = true;
}

/////////////////////////////////////////////////////////////////////////////
// Decomment here to remove circularity criterion from selection procedure //
//lowBoundC = 0;
//highBoundC = 1;
/////////////////////////////////////////////////////////////////////////////

print("Single cell size range = [", lowBound, ", ", highBound, "] pixel^2");
print("Characteristic circularity range = [", lowBoundC, ", ", highBoundC, "]");
if(boundFlag)
	print("WARNING: Circularity bounds have been restored to the average extreme values because of a too low tolerance value");

// Blank images used as 2D arrays: use 32-bit images to have a floating point value for each pixel 
newImage("decoy", "8-bit black", 10, 10, 1); // Robustness with respect to possible empty images
run("Set Measurements...", "area perimeter fit shape feret's redirect=None decimal=5"); // Full set of Descriptors (the same for both soma and nucleus)
run("Measure");
fullheadings = split(String.getResultsHeadings);
headings = newArray();
for (i = 0; i < lengthOf(fullheadings); i++) {
	if (fullheadings[i] != "FeretX" && fullheadings[i] != "FeretY") // Delete uninformative measurements
		headings = Array.concat(headings,fullheadings[i]);
}
selectWindow("decoy");
close();
// MasterMatrix_M to store Cell Shape Descriptors
newImage("MasterMatrix_M", "32-bit white", 1, lengthOf(headings), 1);
// MasterMatrix_N to store Nucleus Shape Descriptors
if(nucleus)
	newImage("MasterMatrix_N", "32-bit white", 1, lengthOf(headings), 1);
// MasterMatrix_O to store Cell and/or Cytoskeleton Orientations
if(orient) {
	newImage("MasterMatrix_O", "32-bit white", 180, 2*sampleNum+3, 1); // +3 because of 1 degree column + 2 separation columns
	selectWindow("MasterMatrix_O");
	for (row = 0; row < 180; row++) { // Add degree column
		setPixel(row,0,-90+row);
	}
}

// Do the Analysis
ScanForAnalysis(input, output, lowBound, highBound, lowBoundC, highBoundC, typicalRadius, headings); // A Morpheus' function - see definition below
print("\nTotal sample files processed: ", shortlist.length); // shortlist.length == sampleNum
selectWindow("MasterMatrix_M");
print("Total cells identified: ", getWidth(), " out of ", objCounter, " detected objects");

// Save MasterMatrix_M as output
flag = saveMatrix("MasterMatrix_M"); // A Morpheus' function - see definition below
// Add custom headings
colname = split(Table.headings);
for (i = 0; i < lengthOf(colname); i++) { // Re-add headings
	Table.renameColumn(colname[i], headings[i]);
}
flagNoUse = File.delete(output + File.separator + "MasterMatrix_M.xls");
saveAs("results", output + File.separator + "MasterMatrix_M.xls");
selectWindow("MasterMatrix_M.xls");
run("Close");

// Do the Nucleus Analysis (Optional)
if(nucleus) {
	N_fileNumber = ScanForNucleus(input, output, suffix, suffix2, identifier, lowBound/10, typicalArea, headings); // A Morpheus' function - see definition below
	print("Total \"", identifier, "\" file processed: ", N_fileNumber);
	selectWindow("MasterMatrix_N");
	if(N_fileNumber > 0)
		print("Total nuclei identified: ", getWidth());
	
	// Save MasterMatrix_N as output
	flag = flag * saveMatrix("MasterMatrix_N"); // Logic AND
	// Add custom headings
	colname = split(Table.headings);
	for (i = 0; i < lengthOf(colname); i++) { // Re-add headings
		Table.renameColumn(colname[i], headings[i]);
	}
	flagNoUse = File.delete(output + File.separator + "MasterMatrix_N.xls");
	saveAs("results", output + File.separator + "MasterMatrix_N.xls");
	selectWindow("MasterMatrix_N.xls");
	run("Close");
}

// Plot and Save MasterMatrix_O as output
if(orient) {
	// Save MasterMatrix_O as output
	selectWindow("MasterMatrix_O");
	run("Flip Horizontally"); // To have positive degrees upward and negative degrees downward
	run("Duplicate...", "title=MasterMatrix_O_image");
	selectWindow("MasterMatrix_O");
	flag = flag * saveMatrix("MasterMatrix_O"); // Logic AND
	// Add custom headings
	colname = split(Table.headings);
	Table.renameColumn(colname[0], "Degree");
	Table.renameColumn(colname[1], "cell->");
	Table.renameColumn(colname[sampleNum+2], "cytoskeleton->");
	for (i = 0; i < sampleNum; i++) {
		Table.renameColumn(colname[i+2], "Cell_"+toString(i+1));
		Table.renameColumn(colname[sampleNum+i+3], "Cyto_"+toString(i+1));
	}
	for (i = 0; i < 180; i++) { // Empty columns (as separators)
		Table.set("cell->", i, "");
		Table.set("cytoskeleton->", i, "");
	}
	Table.showRowNumbers(false);
	flagNoUse = File.delete(output + File.separator + "MasterMatrix_O.xls");
	saveAs("results", output + File.separator + "MasterMatrix_O.xls");
	selectWindow("MasterMatrix_O.xls");
	run("Close");
	
	PlotMM2(sampleNum); // A Morpheus' function - see definition below
}

if(verbose) {
	if(flag == 1) // Logic AND
		print("\nAll output files have been correctly saved to: " + output);
	else
		print("\nWARNING !!! Some error occurred while saving data...");
}

if (isOpen("Results")) {
	selectWindow("Results");
    run("Close");
}

setBatchMode(false); // Note: even in BatchMode (when orient==true) "MasterMatrix_O" is left open because it is the current active image

selectWindow("Log");  // Select Log-window
saveAs("text", output + File.separator + "Log.txt"); // Save session Log

if(orient)
	selectWindow("MasterMatrix_O.tif");

// END of Main

//---------------------------------------------------------------------------------------------------//
// Function Definitions
//---------------------------------------------------------------------------------------------------//

// Search for files with the correct suffix and identifier (e.g. NO DAPI)
// then compute the statistics of interest (i.e. median area and circularity extreme values) from each sample image
function GetSamplingDistributions(input, suffix, suffix2, identifier, lowBound) {

	print("\n");
	list = getFileList(input);
	population = newArray();
	
	for (i = 0; i < list.length; i++) {
		
		// Search for files with .suffix extension and avoid files with "identifier" substring within their name
		if((endsWith(toLowerCase(list[i]), suffix) || endsWith(toLowerCase(list[i]), suffix2)) && indexOf(toLowerCase(list[i]), identifier) == -1) {
			
			open(input + File.separator + list[i]);
			if(verbose)
				print("Pre-Processing: " + input + File.separator + list[i]);
				
			// Segmentation
			MorpheuSegment(list[i], false); // A Morpheus' function - see definition below
			
			// Analyze Particles
			run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel"); // "Click to Remove Scale" option - Force pixel as unit of length
			run("Analyze Particles...", "size=" + lowBound + "-Infinity pixel circularity=0-1 show=Nothing exclude clear include");
			
			objCounter = objCounter + nResults;
			
			// Compute sampling statistics
			if(nResults > 0) {
				shortlist = Array.concat(shortlist, list[i]);
				sampleNum = population.length/3;
				
				// Array of sample areas and sample circularities
				a = newArray();
				c = newArray();
				for (j = 0; j < nResults; j++) {
					a = Array.concat(a, getResult("Area", j));
					c = Array.concat(c, getResult("Circ.", j));
				}
				//Array.print(a); // Debug
				//Array.print(c); // Debug
				
				// Sample statistics of interest
				meda = findMedian(a); // A Morpheus' function - see definition below
				Array.getStatistics(c, minc, maxc);
				stat = newArray(meda, minc, maxc);
				
				if(verbose) {
					print("    Sample_" + sampleNum+1 + " Statistics:");
					print("        Sample Size = ", nResults);
					print("        Median Area = ", d2s(stat[0],1), " pixel^2");
					print("        Circularity Minimum = ", d2s(stat[1],4));
					print("        Circularity Maximum = ", d2s(stat[2],4));
				}
				population = Array.concat(population, stat);
				
				saveAs("tiff", output + File.separator + "Cell" + File.separator + "Img_" + sampleNum+1 + "a - Segmentation - " + list[i]);
				close(); // Close 'Segmentation' image
			}
			else {
				if(verbose)
					print("WARNING: there are no detectable objects in this image. It may be void.");
				close(); // Close 'Segmentation' (empty) image
			}
			
			selectWindow(list[i]);
			close(); // Close original image
		}
	}
	
	return population;
}

// Find Median
function findMedian(x) {
	
	Array.sort(x);
	
	if(x.length == 0)
		median = 0;
	else if ( (x.length % 2) == 1 ) // Odd number of elements
		median = x[(x.length-1)/2];
	else if ( (x.length % 2) == 0 ) // Even number of elements
		median = (x[x.length/2] + x[(x.length/2)-1])/2;
	
	return median;	
}

// Do the Analysis
function ScanForAnalysis(input, output, lowBound, highBound, minCirc, maxCirc, radius, headings) {
	
	if(verbose)
		print("\n");
	
	for (i = 0; i < shortlist.length; i++) {
		if(verbose)
			print("Processing: " + input + File.separator + shortlist[i]);
		
		open(output + File.separator + "Cell" + File.separator + "Img_" + i+1 + "a - Segmentation - " + shortlist[i]);
		
		// Selection of isolated cells
		run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel"); // "Click to Remove Scale" option - Force pixel as unit of length
		run("Analyze Particles...", "size=" + lowBound + "-" + highBound + " pixel circularity=" + minCirc + "-" + maxCirc + " show=Outlines display exclude clear include");
		selectWindow("Img_" + i+1 + "a - Segmentation - " + shortlist[i]);
		run("Analyze Particles...", "size=" + lowBound + "-" + highBound + " pixel circularity=" + minCirc + "-" + maxCirc + " show=Masks display exclude clear include");
		
		if (nResults > 0) {
			selectWindow("Results");
			Table.deleteColumn("FeretX"); // Delete uninformative measurements
			Table.deleteColumn("FeretY"); // Delete uninformative measurements
			saveAs("results", output + File.separator + "Cell" + File.separator + "Descriptors_" + i+1 + ".xls");
		}
		
		selectWindow("MasterMatrix_M");
		if(i == 0)
			extendWidth = nResults;
		else
			extendWidth = getWidth() + nResults;
		
		run("Canvas Size...", "width=" + extendWidth + " height=" + lengthOf(headings) + " position=Center-Left zero");
		for (row = 0; row < nResults; row++) {
			for (col = 0; col < lengthOf(headings); col++)
				setPixel(row+getWidth()-nResults, col, getResult(headings[col], row));
		}
		
		ncells = nResults;
		if(orient) {
			CytoskeletOrient(output, i, shortlist[i], ncells); // A Morpheus' function - see definition below
			if (orientStr == "Double-Scale")
				CellOrient(output, i, radius, shortlist[i], ncells); // A Morpheus' function - see definition below
		}
		
		selectWindow("Drawing of Img_" + i+1 + "a - Segmentation - " + shortlist[i]);
		saveAs("tiff", output + File.separator + "Cell" + File.separator + "Img_" + i+1 + "b - Detection - " + shortlist[i]);
		close(); // Close Cell Detection (Outlines)
		
		selectWindow("Img_" + i+1 + "a - Segmentation - " + shortlist[i]);
		close(); // Close 'Segmentation' image
		
		selectWindow("Mask of Img_" + i+1 + "a - Segmentation - " + shortlist[i]);
		//saveAs("tiff", output + File.separator + "Cell" + File.separator + "Img_" + i+1 + "ab - Mask - " + shortlist[i]); // Debug
		close(); // Close 'Mask' image
	}
	
	//return something;
}

// Do Nucleus Analysis
function ScanForNucleus(input, output, suffix, suffix2, identifier, lowBound, highBound, headings) {
	
	list = getFileList(input);
	j = 0; // Count just the files with the correct suffix and identifier
	
	if(verbose)
		print("\n");
	
	for (i = 0; i < list.length; i++) {
		if((endsWith(toLowerCase(list[i]), suffix) || endsWith(toLowerCase(list[i]), suffix2)) && indexOf(toLowerCase(list[i]), identifier) != -1) { // Search for files with the "identifier" substring within their name
			
			open(input + File.separator + list[i]);
			if(verbose)
				print("Processing: " + input + File.separator + list[i]);
			
			// Segmentation
			MorpheuSegment(list[i], true); // A Morpheus' function - see definition below

			// Analyze Particles
			run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel"); // "Click to Remove Scale" option - Force pixel as unit of length
			run("Analyze Particles...", "size=" + lowBound + "-" + highBound + " pixel circularity=0.7-1 show=Outlines display exclude clear include");

			if (nResults > 0) {
				selectWindow("Results");
				Table.deleteColumn("FeretX"); // Delete uninformative measurements
				Table.deleteColumn("FeretY"); // Delete uninformative measurements
				saveAs("results", output + File.separator + "Nucleus" + File.separator + "N_Descriptors_" + j+1 + ".xls");
			}
			
			selectWindow("MasterMatrix_N");
			if(j == 0)
				extendWidth = nResults;
			else
				extendWidth = getWidth() + nResults;
			
			run("Canvas Size...", "width=" + extendWidth + " height=" + lengthOf(headings) + " position=Center-Left zero");
			for (row = 0; row < nResults; row++) {
				for (col = 0; col < lengthOf(headings); col++)
					setPixel(row+getWidth()-nResults, col, getResult(headings[col], row));
			}
			
			selectWindow("Drawing of Segment");
			saveAs("tiff", output + File.separator + "Nucleus" + File.separator + "N_Img_" + j+1 + "b - Detection - " + list[i]);
			close(); // Close Nuclei Detection (Outlines)
			
			selectWindow("Segment");
			saveAs("tiff", output + File.separator + "Nucleus" + File.separator + "N_Img_" + j+1 + "a - Segmentation - " + list[i]);
			close(); // Close Segmentation image
			
			selectWindow(list[i]);
			close(); // Close original image
			
			j = j+1;
		}
	}
	
	return j;
}

// Segmentation Procedure
function MorpheuSegment(inName, isNucleus) {
	
	selectWindow(inName);
	run("Duplicate...", "title=Segment");
	
	// Enhance contrast to increase the dynamic range
	MorpheusEnhance(); // A Morpheus' function - see definition below
	
	rollingRad = ((getWidth()+getHeight())/2)/6; // Rolling ball diameter == 1/3 of image linear dimension
	run("Subtract Background...", "rolling=" + rollingRad + " sliding");
	run("Smooth");
	run("Auto Threshold", "method=Li");
	run("Make Binary");
	if (isNucleus)
		run("Watershed");
	
	//return something;
}

// Robust Enhance Contrast
function MorpheusEnhance() {
	
	run("8-bit"); // RGB/16-bit to 8-bit conversion
	
	minLev = 64; // Minimum number of gray levels required for saturation
	sat = 0.05; // Starting percentage of saturated pixels
	getStatistics(size, meanNoUse, minPix, maxPix); // size == getWidth()*getHeight()
	getHistogram(values, counts, 256);
	
	if ((maxPix-minPix) > minLev) {
		do {
			cum = 0; // "Right" cumulative distribution function (1-CDF)
			i = 0; // Threshold "right" quantile
			
			do {
				cum = cum + counts[255-i];
				i = i + 1;
			} while (cum < (sat/100)*size);
			
			// Decrease sat when gray levels are too few
			//(it avoids background-noise amplification by false contours)
			sat = sat / 2;
			
		} while (256-i < minLev);
	} else {
		i = 256 - maxPix;
		sat = 0;
	}
	
	maxLocs = Array.findMaxima(counts, sqrt(size)); // Returns an array holding the peak positions sorted with descending strength
	mode = maxLocs[0];
	
	// Debug
	/*
	print("Max set to: " + 256-i);
	print("Up-Saturated: " + sat*2 + " %");
	print("Histogram peaks:");
	Array.print(maxLocs);
	print("Mode: " + mode);
	print("Min set to: " + mode+1);
	*/
	
	setMinAndMax(mode+1, 256-i); // min = mode+1 to exclude background
	run("Apply LUT");
	
	//return something;
}

// Cytoskeleton Orientation
function CytoskeletOrient(output, serial, file, ncells) {
	
	open(input + File.separator + file);
	
	run("Duplicate...", "title=Original");
	run("Sharpen"); // Optional
	
	// Enhance contrast to increase the dynamic range
	MorpheusEnhance(); // A Morpheus' function - see definition below
	
	if (orientStr == "All Filaments") {
		selectWindow("Img_" + i+1 + "a - Segmentation - " + file); // Morpheus Segmentation Mask
		run("Duplicate...", "title=MultiplMask");
		run("Invert");
	} else { //if (orientStr == "Double-Scale")
		selectWindow("Mask of Img_" + i+1 + "a - Segmentation - " + file); // Morpheus Detection Mask
		run("Duplicate...", "title=MultiplMask");
	}
	
	run("Divide...", "value=255"); // {0,1} Values (Multiplicative Mask)
	imageCalculator("Multiply create", "Original", "MultiplMask"); // This creates a new image called "Result of Original"
	
	// gradient=4 -> Gaussian Gradient
	run("OrientationJ Distribution", "tensor=1.0 gradient=4 orientation=on coherency=on radian=off binary=on histogram=on table=on min-coherency=1.0 min-energy=1.0");
	
	// Normalize Orientation Histogram and fill out MasterMatrix_O
	selectWindow("OJ-Distribution-1");
	orVal = Table.getColumn("Slice1"); // Array
	run("Close"); // Close Table
	
	selectWindow("MasterMatrix_O");
	sampleNum = (getHeight()-3)/2;
	
	if (orientStr == "All Filaments") {
		for (row = 0; row < 180; row++) {
			setPixel(row, sampleNum+3+serial, orVal[row]); // No Normalization
		}
	} else { //if (orientStr == "Double-Scale")
		norm = 0;
		for (row = 0; row < 180; row++) {
			norm = norm + orVal[row];
		}
		for (row = 0; row < 180; row++) {
			setPixel(row, sampleNum+3+serial, (orVal[row]/norm)*ncells);
		}
	}
	
	// Save Histogram
	selectWindow("OJ-Histogram-1-slice-1");
	if(full) {
		saveAs("tiff", output + File.separator + "Cell" + File.separator + "Img_" + serial+1 + "g - CytoHisto - " + file);
		close();
	}
	while(isOpen("OJ-Histogram-1-slice-1")) { // To close multiple Histogram windows
		selectWindow("OJ-Histogram-1-slice-1");
		run("Close");
	}
	
	// Build the Color Survey
	makeColorSurvey("Result of Original"); // A Morpheus' function - see definition below
	
	selectWindow("OJ-Orientation-1");
	run("Close"); // Don't know why, but OJ output images need to be closed by run("Close") when running in BatchMode. close() command is not effective...
	selectWindow("OJ-Coherency-1");
	run("Close");
	
	selectWindow(file);
	close(); // Close original
	selectWindow("Original");
	close(); // Close original duplicate
	selectWindow("MultiplMask");
	close(); // Close multiplicative Mask
	selectWindow("Result of Original");
	close(); // Close the result of imageCalculator("Multiply create")
	
	selectWindow("ColorSurvey");
	saveAs("tiff", output + File.separator + "Cell" + File.separator + "Img_" + serial+1 + "d - CytoOrient - " + file);
	close(); // Close color survey
	
	selectWindow("OJ-Orientation Mask-1");
	if(full) { // Useful for high-energy noise artifacts and edge effect checking
		
		run("16-bit"); // Dynamic range -> [0,65535]
		run("Spectrum");
		run("Apply LUT");
		run("RGB Color"); // Dynamic range -> [0,255]
		run("HSB Stack"); // Dynamic range -> [0,255]
		selectWindow("OJ-Binary Mask-1");
		run("Select All");
		run("Copy");
		selectWindow("OJ-Orientation Mask-1");
		setSlice(3); // Brightness
		run("Paste");
		run("RGB Color");
		
		saveAs("tiff", output + File.separator + "Cell" + File.separator + "Img_" + serial+1 + "e - S-CytoOrient - " + file);
	}
	run("Close"); // OJ-Orientation Mask
	
	selectWindow("OJ-Binary Mask-1");
	run("Close");
	
	//return something;
}

// Cell Orientation
function CellOrient(output, serial, radius, file, ncells) {
	
	selectWindow("Mask of Img_" + i+1 + "a - Segmentation - " + file);
	// gradient=4 -> Gaussian Gradient
	run("OrientationJ Distribution", "tensor=" + radius + " gradient=4 orientation=on coherency=on radian=off histogram=on table=on min-coherency=0.0 min-energy=0.0 ");
	
	// Close table and histograms
	if (isOpen("OJ-Distribution-1")) {
		selectWindow("OJ-Distribution-1");
		run("Close");
	}
	while(isOpen("OJ-Histogram-1-slice-1")){ // To close multiple Histogram windows
		selectWindow("OJ-Histogram-1-slice-1");
		run("Close");
	}
	
	// Build the Color Survey
	makeColorSurvey("Mask of Img_" + i+1 + "a - Segmentation - " + file); // A Morpheus' function - see definition below
	
	selectWindow("ColorSurvey");
	saveAs("tiff", output + File.separator + "Cell" + File.separator + "Img_" + serial+1 + "c - CellOrient - " + file);
	close(); // Close saved tiff file
	
	// Build Coherency-weighted orientation histogram
	selectWindow("Mask of Img_" + i+1 + "a - Segmentation - " + file);
	run("Duplicate...", "title=dup-for-bright"); // Use Brightness-slice -> Original Binarized Image as Mask
	run("8-bit"); // Safety
	run("Divide...", "value=255"); // {0,1} Values (Multiplicative Mask)
	
	setOption("ScaleConversions", false); // Don't scale when converting to 8-bit
	
	selectWindow("OJ-Coherency-1");
	run("Multiply...", "value=255"); // [0,1] -> [0,255]
	run("8-bit");
	imageCalculator("Multiply create", "OJ-Coherency-1", "dup-for-bright"); // Mask Coherency to make weights for histogram
	selectWindow("Result of OJ-Coherency-1");
	rename("weight");
	
	selectWindow("dup-for-bright");
	close();
	
	selectWindow("OJ-Orientation-1");
	run("Duplicate...", "title=dup-for-hue2");
	run("Add...", "value=90"); // [-90,+90] -> [0,180]
	run("8-bit");
	//run("Add Specified Noise...", "standard=1"); // Noise will smooth histogram to avoid quantization spikes
	
	setOption("ScaleConversions", true); // Restore default
	
	selectWindow("OJ-Orientation-1");
	run("Close");
	selectWindow("OJ-Orientation Mask-1");
	run("Close");
	selectWindow("OJ-Coherency-1");
	run("Close");
	
	selectWindow("weight");
	getStatistics(areaW, meanW, minW, maxW);
	step = floor(maxW/10); // Just 10-level weights (10 instead of 256 to speed up the computation)
	//print("Max of Coherency = " + maxW); // Debug
	//print("Ten-level step = " + step); // Debug
	
	mainhisto = newArray(180);
	
	if(step >= 1) { // ...that is (maxW >= 10)
		for(h = 0; h < 10; h++) {
			
			// Debug
			/*
			if(h+1 == 10)
				print("level", h+1, ":", h*step+1, "-", maxW);
			else
				print("level", h+1, ":", h*step+1, "-", (h+1)*step, "\n");
			*/
			
			selectWindow("weight");
			run("Duplicate...", "title=temp");		
			setAutoThreshold("Default dark");
			call("ij.plugin.frame.ThresholdAdjuster.setMode", "B&W");
			
			if(h+1 == 10)
				setThreshold(h*step+1, maxW);
			else
				setThreshold(h*step+1, (h+1)*step);
			
			run("Convert to Mask");	
			run("Divide...", "value=255"); // {0,1} Values (Multiplicative Mask)
			
			imageCalculator("Multiply create", "temp", "dup-for-hue2");
			selectWindow("Result of temp");
			getStatistics(areaNoUse, meanNoUse, minNoUse, maxNoUse, stdNoUse, histo);
			//getHistogram(0, counts, 256); // As an alternative, to return just the histogram
			
			selectWindow("Result of temp");
			close();
			selectWindow("temp");
			close();
			
			histo[0] = 0; // Delete background bin (and zero-Coherency pixels)
			for(k = 0; k < 180; k++){
				mainhisto[k] = mainhisto[k] + histo[k]*((h+1)*step-step/2); // Build weighted histogram
			}
		}
	}
	else{
		//print("unique level", ":", 1, "-", maxW); // Debug
		
		selectWindow("Mask of Img_" + i+1 + "a - Segmentation - " + file);
		run("Duplicate...", "title=temp");
		run("Divide...", "value=255"); // {0,1} Values (Multiplicative Mask)
		imageCalculator("Multiply create", "temp", "dup-for-hue2");
		selectWindow("Result of temp");
		getStatistics(areaNoUse, meanNoUse, minNoUse, maxNoUse, stdNoUse, histo);
		
		selectWindow("Result of temp");
		close();
		selectWindow("temp");
		close();
		
		histo[0] = 0; // Delete background bin (and zero-Coherency pixels)
		for(k = 0; k < 180; k++){
			mainhisto[k] = histo[k]*5; // Build weighted histogram (5 is the midpoint in the coherency range [1-9])
		}
	}
	
	selectWindow("dup-for-hue2");
	close();
	selectWindow("weight");
	close();
	
	mainhisto[0] = round((mainhisto[1] + mainhisto[179])/2); // Circular junction
	
	if(full) {
		Plot.create("MainHistogram", "Orientation in Degrees", "Distribution of orientation", degreeAxis, mainhisto); // Same labels of OrientationJ
		Plot.show();
		selectWindow("MainHistogram");
		saveAs("tiff", output + File.separator + "Cell" + File.separator + "Img_" + serial+1 + "f - CellHisto - " + file);
		close();
	}
	
	// Normalization factor
	norm = 0;
	for (row = 0; row < 180; row++) {
		norm = norm + mainhisto[row];
	}
	selectWindow("MasterMatrix_O");
	for (row = 0; row < 180; row++) {
		setPixel(row, serial+2, (mainhisto[row]/norm)*ncells);
	}
	
	//return something;
}

// Plot MasterMatrix_O
function PlotMM2(sampleNum) {
	
	// Transpose Matrix
	selectWindow("MasterMatrix_O_image");
	run("Rotate 90 Degrees Right");
	run("Flip Horizontally");
	
	makeRectangle(2, 0, 2*sampleNum+1, 180); // Crop "Degree" column
	run("Crop");
	
	getStatistics(areaM, meanM, minM, maxM);
	setMinAndMax(0, maxM);
	run("16-bit"); // Dynamic range -> [0,65535]
	
	makeRectangle(0, 0, sampleNum, 180); // Enhance Contrast to normalize Soma and Cytoskeleton separately, allowing for comparison
	run("Enhance Contrast", "saturated=0");
	run("Apply LUT"); // Dynamic range -> [0,65535]
	makeRectangle(sampleNum+1, 0, sampleNum, 180);
	run("Enhance Contrast", "saturated=0");
	run("Apply LUT"); // Dynamic range -> [0,65535]
	run("Select None");
	
	// Draw a separation white line between Shape and Cytoskeleton Orientations and apply a LUT
	getStatistics(areaMM, meanMM, minMM, maxMM);
	for (row = 0; row < 180; row++) {
		setPixel(sampleNum, row, maxMM);
	}
	run("Fire");
	
	newHeight = 4*getHeight(); // Magnification 4x
	newWidth = 4*getWidth();
	run("Size...", "width=" + newWidth + " height=" + newHeight + " constrain average interpolation=None");
	run("Canvas Size...", "width=" + getWidth()+60 + " position=Center-Right zero");
	
	// Draw degree vertical axis
	newImage("VerticalAxis", "16-bit black", 60, 720, 1);
	setColor("white");
	setFont("SansSerif", 16, "bold");
	setJustification("right");
	setLineWidth(2);
	drawLine(55, 0, 55, 720); // Vertical Axis
	for (k = 1; k <= 18; k++) {
		drawLine(55, k*4*10-1-2, 50, k*4*10-1-2); // Horizontal Ticks
		if(k < 9)
			drawString("+ " + abs(10*(9-k)), 45, k*4*10-1-2); // Positive labels (-1 to account for circularity constraint; -2 to center with respect to each "4x4-pixel")
		else if (k > 9)
			drawString("- " + abs(10*(9-k)), 45, k*4*10-1-2); // Negative labels
		else
			drawString(abs(10*(9-k)), 45, k*4*10-1-2); // Zero label
	}
	run("Invert");
	
	makeRectangle(0, 0, 60, 720);
	run("Copy");
	selectWindow("MasterMatrix_O_image");
	makeRectangle(0, 0, 60, 720);
	run("Paste");
	selectWindow("VerticalAxis");
	close();
	
	selectWindow("MasterMatrix_O_image");
	run("Select None");
	saveAs("tiff", output + File.separator + "MasterMatrix_O.tif");
	
	//return something;
}

// Print General Log Information
function systeminfo(morphversion) {
	
	print("\\Clear"); // Empty the Log
	print(">------------------------<");
	print("  Starting Morpheus");
	print(">------------------------<");
	
	// Print the Date
	MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
	DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	TimeString = "Date: " + DayNames[dayOfWeek] + " ";
	if (dayOfMonth<10) {TimeString = TimeString + "0";}
	TimeString = TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+"\nTime: ";
	if (hour<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+hour+":";
	if (minute<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+minute+":";
	if (second<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+second;
	print("\n" + TimeString);
	
	// Print Versioning Info
	print("\nArchitecture:   " + getInfo("os.arch"));
	print("Operating system:   " + getInfo("os.name") + " (version " + getInfo("os.version") + ")");
	print("Java version:   " + getInfo("java.version"));
	print(getInfo("java.runtime.name"));
	print("Java Virtual Machine:   " + getInfo("java.vm.name") + " (build " + getInfo("java.vm.version") + ")");
	print("ImageJ version:   " + getVersion());
	
	print("\nMorpheus Version:   " + morphversion);
	print("    Anti-spot lower bound = " + lowBound + " pixel^2");
	print("    Tolerance = " + tol + " sigma");
	if(orient) {print("    Orientation analysis: " + orientStr + " mode");}
	if(nucleus) {print("    Nucleus analysis: identifier = " + identifier + ")");}
	if(full) {print("    Full output option");}
	if(verbose) {print("    Verbose log option");}
	
	print("\nRequired Plugin Versions: ");
	pluginlist = getFileList(getDirectory("plugins"));
	for (i = 0; i < pluginlist.length; i++) {
		if(endsWith(pluginlist[i], ".jar")) {
			if((indexOf(toLowerCase(pluginlist[i]),"auto_threshold") != -1) || (indexOf(toLowerCase(pluginlist[i]),"orientationj") != -1))
				print("    " + pluginlist[i]);
		}
	}
	
	//return something;
}

// Build the Color Survey in HSB space and then convert to RGB
function makeColorSurvey(brightnessName) {
	
	selectWindow("OJ-Orientation-1");
	newImage("ColorSurvey", "RGB black", getWidth(), getHeight(), 1);
	run("HSB Stack");
	
	setOption("ScaleConversions", false); // Don't scale when converting to 8-bit
	
	selectWindow("OJ-Orientation-1");
	run("Duplicate...", "title=dup-for-hue"); // Rescale [-90,+90] -> [0,255]
		run("Add...", "value=90");
		run("Divide...", "value=180");
		run("Multiply...", "value=255");
		run("8-bit");
	run("Select All");
	run("Copy");
	selectWindow("ColorSurvey");
	setSlice(1); // Hue
	run("Paste");
	
	selectWindow("OJ-Coherency-1");
	run("Duplicate...", "title=dup-for-sat"); // Rescale [0,1] -> [0,255]
		run("Multiply...", "value=255");
		run("8-bit");
	run("Select All");
	run("Copy");
	selectWindow("ColorSurvey");
	setSlice(2); // Saturation
	run("Paste");
	
	setOption("ScaleConversions", true); // Restore default
	
	selectWindow(brightnessName);
	run("Select All");
	run("Copy");
	selectWindow("ColorSurvey");
	setSlice(3); // Brightness
	run("Paste");
	
	run("RGB Color");
	
	selectWindow("dup-for-hue");
	close();
	selectWindow("dup-for-sat");
	close();
	
	//return something;
}

function saveMatrix(inName) {
	
	// Transpose Matrix
	selectWindow(inName);
	run("Rotate 90 Degrees Right");
	run("Flip Horizontally");
	
	// Save
	saveAs("text image", output + File.separator + inName + ".txt"); // ...can't save in .xls directly
	close();
	flagNoUse = File.delete(output + File.separator + inName + ".xls"); // Delete any files with the same name, produced by possible previous runs
	flag = File.rename(output + File.separator + inName + ".txt", output + File.separator + inName + ".xls"); // Change extension - Returns "1" (true) if successful
	Table.open(output + File.separator + inName + ".xls");
	
	return flag;
}

// Morpheus default options
function setDefaultOptions() {
	
	setOption("BlackBackground", false); // Binary and Threshold option
	//setOption("InvertY", false); // It requires an open image
	setOption("JFileChooser", false);
	setOption("Show All", true);
	run("Input/Output...", "jpeg=85 gif=-1 file=.csv use_file copy_row save_column save_row");
		setOption("ShowRowNumbers", true); // For "Results" window
		Table.showRowNumbers(true); // For tables
	setOption("ScaleConversions", true); // Scale from min–max to 0–255 when converting from 32– or 16-bit to 8–bit (Morpheus' default)
	//run("Profile Plot Options...", "width=350 height=200 draw"); // To fix plot size (always 350×200 pixels) and style
	
	//return something;
}
