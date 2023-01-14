// MORPHEUS -- release 2018
// ImageJ 1.x macro language version

// Universal I/O parameter notation
#@ File(label="Input directory", style="directory") input
#@ File(label="Output directory", style="directory") output
#@ String(label="File suffix", value=".tif") suffix
#@ String(label="Nucleus file identifier", value="dapi") identifier
#@ Integer(label="Anti-Spot Lower Bound (pixel^2)", value=200) lowBound
#@ String(label="Tolerance (sigma)", choices={"1","2","3","4","5","6"}, style="listBox", value="3") tol
#@ Boolean(label="Orientation Analysis", value=false) orient
#@ Boolean(label="Nucleus Analysis", value=false) nucleus
#@ Boolean(label="Full Output", value=false) full
#@ Boolean(label="Verbose Log", value=true) verbose

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
if(lengthOf(identifier) == 0) {
	identifier = "dapi";
}

// Organize Output Directory
File.makeDirectory(output + File.separator + "Cell");
if(nucleus) {
	File.makeDirectory(output + File.separator + "Nucleus");
}

// Tolerance -> string to integer
tol = parseInt(tol);

// Starting Morpheus Macro
print("\n");
print(">------------------------<");
print("  Starting Morpheus");
print(">------------------------<");
print("\n");

// Batch Mode! ...to speed up the macro
setBatchMode(true);

// Global variables
var objCounter = 0; //Count all the objects detected

// Learn the characteristic features of the single cell from the whole dataset
run("Set Measurements...", "area shape redirect=None decimal=3");
// "pop" is an array featuring median areas and extreme circularities of each sample image alternately
pop = GetSamplingDistributions(input, suffix, suffix2, identifier, lowBound); // A function of mine - see below

// Check for contents
if(pop.length == 0) {
	messageOut = "No file with " + suffix + " suffix and suitable for cell morphometry is present in the folder\n" + input;
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

// Define typical cell area and circularity range as the means of their corresponding sampling distributions 
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
//if(0.5*(typicalArea-tol*AreaStd/sqrt(sampleNum)) > lowBound) {
//	lowBound = 0.5*(typicalArea-tol*AreaStd/sqrt(sampleNum));
//}
// Highest tolerance (6):	typicalMinC + 0*MinCStd AND typicalMaxC - 0*MaxCStd
// Mid tolerance (3):		typicalMinC + 3*MinCStd AND typicalMaxC - 3*MaxCStd
// Lowest tolerance (1):	typicalMinC + 5*MinCStd AND typicalMaxC - 5*MaxCStd
lowBoundC = typicalMinC + (6-tol)*MinCStd/sqrt(sampleNum);
highBoundC = typicalMaxC - (6-tol)*MaxCStd/sqrt(sampleNum);
boundFlag = false;
if(highBoundC <= lowBoundC) {
	lowBoundC = typicalMinC;
	highBoundC = typicalMaxC;
	boundFlag = true;
}

print("Single cell size range = [", lowBound, ", ", highBound, "] pixel^2");
print("Characteristic circularity range = [", lowBoundC, ", ", highBoundC, "]");
if(boundFlag) {
	print("WARNING: Circularity bounds have been restored to the average extreme values because of a too low tolerance value");
}

// Blank images used as 2D arrays: use 32-bit images to have a floating point value for each pixel 
// MegaMatrix1 to store Shape Descriptors
run("Summarize"); // Display Results
run("Set Measurements...", "area perimeter bounding fit shape feret's redirect=None decimal=3"); // Full set of Descriptors (the same both for soma and nucleus)
headings = split(String.getResultsHeadings);
newImage("MegaMatrix1", "32-bit white", 1, lengthOf(headings), 1);
// MegaMatrix1_N to store Nucleus Shape Descriptors
if(nucleus) {
	newImage("MegaMatrix1_N", "32-bit white", 1, lengthOf(headings), 1);
}
// MegaMatrix2 to store Shape and Cytoskeleton Orientations
if(orient) {
	newImage("MegaMatrix2", "32-bit white", 180, 2*sampleNum+3, 1);
	selectWindow("MegaMatrix2");
	for (row = 0; row < 180; row++) {
		setPixel(row,0,-90+row);
	}
}

// Do the Analysis
fileNumber = ScanForAnalysis(input, output, suffix, suffix2, identifier, lowBound, highBound, lowBoundC, highBoundC, typicalRadius, headings); // A function of mine - see below
print("\nTotal sample files processed: ", fileNumber);
selectWindow("MegaMatrix1");
print("Total cells identified: ", getWidth(), " out of ", objCounter, " detected objects");

run("Rotate 90 Degrees Right"); // Transpose MegaMatrix1
run("Flip Horizontally");
saveAs("text image", output + File.separator + "MegaMatrix1.txt"); // ...can't save in .xls directly
close();
flagNoUse = File.delete(output + File.separator + "MegaMatrix1.xls"); // Delete any files with the same name, produced by possible previous runs
flag = File.rename(output + File.separator + "MegaMatrix1.txt", output + File.separator + "MegaMatrix1.xls"); // Returns "1" (true) if successful

if(nucleus) {
	N_fileNumber = ScanForNucleus(input, output, suffix, suffix2, identifier, typicalArea/10, typicalArea, headings); // A function of mine - see below - WARNING!! C'è accordo sui bounds?
	print("\nTotal \"", identifier, "\" file processed: ", N_fileNumber);
	selectWindow("MegaMatrix1_N");
	if(N_fileNumber > 0) {
		print("Total nuclei identified: ", getWidth());
	}
	
	run("Rotate 90 Degrees Right"); // Transpose MegaMatrix1_N
	run("Flip Horizontally");
	saveAs("text image", output + File.separator + "MegaMatrix1_N.txt"); // ...can't save in .xls directly
	close();
	flagNoUse = File.delete(output + File.separator + "MegaMatrix1_N.xls"); // Delete any files with the same name, produced by possible previous runs
	flag = flag * File.rename(output + File.separator + "MegaMatrix1_N.txt", output + File.separator + "MegaMatrix1_N.xls"); // Logic AND
}

if(orient) {
	selectWindow("MegaMatrix2");
	run("Rotate 90 Degrees Right"); // Transpose MegaMatrix2
	run("Flip Horizontally");
	run("Flip Vertically"); // To have positive degrees upward and negative degrees downward
	saveAs("text image", output + File.separator + "MegaMatrix2.txt"); // ...can't save in .xls directly
	flagNoUse = File.delete(output + File.separator + "MegaMatrix2.xls"); // Delete any files with the same name, produced by possible previous runs
	flag = flag * File.rename(output + File.separator + "MegaMatrix2.txt", output + File.separator + "MegaMatrix2.xls"); // Logic AND

	PlotMM2(sampleNum); // A function of mine - see below	
}

if(verbose) {
	if(flag == 1)
		print("\nAll output files have been correctly saved to: " + output);
	else
		print("\nWARNING !!! Some error occurred while saving data...");
}

setBatchMode(false); // Note: if orient==true "MegaMatrix2" will not be closed because it is the active image

//---------------------------------------------------------------------------------------------------//
// Function Definitions
//---------------------------------------------------------------------------------------------------//

// Scan to find files with the correct suffix and identifier (e.g. NO DAPI)
// then compute the statistics of interest (i.e. median area and circularity extreme values) from each sample image
function GetSamplingDistributions(input, suffix, suffix2, identifier, lowBound) {
	
	list = getFileList(input);
	
	population = newArray();
	
	for (i = 0; i < list.length; i++) {
		// Search for files with .suffix extension and avoid files with "identifier" substring within their name
		if((endsWith(toLowerCase(list[i]), suffix) || endsWith(toLowerCase(list[i]), suffix2)) && indexOf(toLowerCase(list[i]), identifier) == -1) {
			
			open(input + File.separator + list[i]);
			if(verbose) {
				print("Pre-Processing: " + input + File.separator + list[i]);
			}
			
			run("Duplicate...", "title=duplicate.tif");
			run("8-bit");
			run("Auto Threshold", "method=Li BlackBackground=false");
			run("Make Binary");
			run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel"); // "Click to Remove Scale" option - Force pixel as unit of length
			run("Analyze Particles...", "size=" + lowBound + "-Infinity pixel show=Nothing exclude clear include");
			
			objCounter = objCounter + nResults;
			
			// Array of sample areas
			a = newArray();
			for (j = 0; j < nResults; j++) {
				a = Array.concat(a, getResult("Area", j));
			}
			//Array.print(a); // Just for testing purpose
			
			// Array of sample circularities
			c = newArray();
			for (j = 0; j < nResults; j++) {
				c = Array.concat(c, getResult("Circ.", j));
			}
			//Array.print(c); // Just for testing purpose
			
			// Statistics of interest
			meda = findMedian(a); // A function of mine - see below
			Array.getStatistics(c, minc, maxc);
			stat = newArray(meda, minc, maxc);
			if(verbose) {
				print("Sample_" + (population.length/3)+1 + " Statistics:");
				print("   Sample Size = ", nResults);
				print("   Median Area = ", stat[0], " pixel^2");
				print("   Circularity Minimum = ", stat[1]);
				print("   Circularity Maximum = ", stat[2]);
			}
			
			selectWindow("duplicate.tif");
			close();
			selectWindow(list[i]);
			close();
			
			population = Array.concat(population, stat);
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
function ScanForAnalysis(input, output, suffix, suffix2, identifier, lowBound, highBound, minCirc, maxCirc, radius, headings) {
	
	list = getFileList(input);
	j = 0; // Count only the files with the correct suffix and identifier
	
	if(verbose) {
		print("\n");
	}
	
	for (i = 0; i < list.length; i++) {
		if((endsWith(toLowerCase(list[i]), suffix) || endsWith(toLowerCase(list[i]), suffix2)) && indexOf(toLowerCase(list[i]), identifier) == -1) { // Avoid files with the "identifier" substring within their name
			
			open(input + File.separator + list[i]);
			if(verbose) {
				print("Processing: " + input + File.separator + list[i]);
			}
			
			run("Duplicate...", "title=morfo.tif");
			
			if(orient) {
				CytoskeletOrient(output, j, list[i]); // A function of mine - see below
			}
			
			selectWindow("morfo.tif");
			run("8-bit");
			run("Auto Threshold", "method=Li BlackBackground=false");
			run("Make Binary");
			run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel"); // "Click to Remove Scale" option - Force pixel as unit of length
			run("Analyze Particles...", "size=" + lowBound + "-" + highBound + " pixel circularity=" + minCirc + "-" + maxCirc + " show=Outlines exclude clear include");
			
			saveAs("results", output + File.separator + "Cell" + File.separator + "Descriptors_" + j+1 + ".xls");
			
			selectWindow("MegaMatrix1");
			
			if(j == 0)
				extendWidth = nResults;
			else
				extendWidth = getWidth() + nResults;
			
			run("Canvas Size...", "width=" + extendWidth + " height=" + lengthOf(headings) + " position=Center-Left zero");	
			for (row = 0; row < nResults; row++) {
				for (col = 0; col < lengthOf(headings); col++)
					setPixel(row+getWidth()-nResults, col, getResult(headings[col], row));
			}
			
			if(orient) {
				CellOrient(output, j, radius, list[i]); // A function of mine - see below
			}
			
			selectWindow("Drawing of morfo.tif");
			saveAs("tiff", output + File.separator + "Cell" + File.separator + "Img_" + j+1 + "b - Outlines - " + list[i]);
			close(); // Close Cell Outlines
			
			selectWindow("morfo.tif");
			saveAs("tiff", output + File.separator + "Cell" + File.separator + "Img_" + j+1 + "a - Thresholded - " + list[i]);
			close(); // Close Thresholded image
			
			selectWindow(list[i]);
			close(); // Close original image
			
			j = j+1;
		}
	}
	
	return j;
}

// Do Nucleus Analysis
function ScanForNucleus(input, output, suffix, suffix2, identifier, lowBound, highBound, headings) {
	
	list = getFileList(input);
	j = 0; // Count only the files with the correct suffix and identifier
	
	if(verbose) {
		print("\n");
	}
	
	for (i = 0; i < list.length; i++) {
		if((endsWith(toLowerCase(list[i]), suffix) || endsWith(toLowerCase(list[i]), suffix2)) && indexOf(toLowerCase(list[i]), identifier) != -1) { // Search for files with the "identifier" substring within their name
			
			open(input + File.separator + list[i]);
			if(verbose) {
				print("Processing: " + input + File.separator + list[i]);
			}
			
			run("Duplicate...", "title=morfo.tif");
			run("8-bit");
			run("Auto Threshold", "method=Li BlackBackground=false");
			run("Make Binary");
			run("Watershed");
			run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel"); // "Click to Remove Scale" option - Force pixel as unit of length
			run("Analyze Particles...", "size=" + lowBound + "-" + highBound + " pixel show=Outlines exclude clear include");
			
			saveAs("results", output + File.separator + "Nucleus" + File.separator + "N_Descriptors_" + j+1 + ".xls");
			
			selectWindow("MegaMatrix1_N");
			
			if(j == 0)
				extendWidth = nResults;
			else
				extendWidth = getWidth() + nResults;
			
			run("Canvas Size...", "width=" + extendWidth + " height=" + lengthOf(headings) + " position=Center-Left zero");
			for (row = 0; row < nResults; row++) {
				for (col = 0; col < lengthOf(headings); col++)
					setPixel(row+getWidth()-nResults, col, getResult(headings[col], row));
			}
			
			selectWindow("Drawing of morfo.tif");
			saveAs("tiff", output + File.separator + "Nucleus" + File.separator + "N_Img_" + j+1 + "b - Outlines - " + list[i]);
			close(); // Close Nuclei Outlines
			
			selectWindow("morfo.tif");
			saveAs("tiff", output + File.separator + "Nucleus" + File.separator + "N_Img_" + j+1 + "a - Thresholded - " + list[i]);
			close(); // Close Thresholded image
			
			selectWindow(list[i]);
			close(); // Close original image
			
			j = j+1;
		}
	}
	
	return j;
}

//Cytoskeleton Orientation
function CytoskeletOrient(output, serial, file) {
	
	selectWindow("morfo.tif");
	run("Duplicate...", "title=morfo_enhance.tif");
	run("Enhance Contrast", "saturated=0.35"); // Auto Brightness/Contrast
	run("Apply LUT");
	run("OrientationJ Distribution", "log=0.0 tensor=1.0 gradient=0 min-coherency=0.0 min-energy=1.0 harris-index=on color-survey=on s-color-survey=on s-distribution=on hue=Orientation sat=Coherency bri=Original-Image");
	
	selectWindow("morfo_enhance.tif");
	close(); // Close duplicate
	
	selectWindow("S-Distribution-1");
	Plot.showValues();
	//saveAs("results", output + File.separator + "Cell" + File.separator + "Cytoskeleton_" + serial+1 + ".xls"); // Decomment this line to save a single Excel file for every sample image
	close(); // Close S-Distribution (histogram) window
	
	selectWindow("MegaMatrix2");
	sampleNum = (getHeight()-3)/2;
	for (row = 0; row < nResults; row++) {
		setPixel(row, sampleNum+3+serial, getResult("Y", row));
	}
	
	selectWindow("Color-survey-1");
	saveAs("tiff", output + File.separator + "Cell" + File.separator + "Img_" + serial+1 + "d - Cytoskeleton - " + file);
	close(); // Close color survey
	
	selectWindow("S-Color-survey-1");
	if(full) { // Useful for edge effect checking
		saveAs("tiff", output + File.separator + "Cell" + File.separator + "Img_" + serial+1 + "e - S-Cytoskeleton - " + file);
	}
	close(); // Close Selection color survey
	
	//return something;
}

// Cell Orientation
function CellOrient(output, serial, radius, file) {
	
	selectWindow("morfo.tif");
	run("OrientationJ Distribution", "log=0.0 tensor=" + radius + " gradient=0 min-coherency=0.0 min-energy=0.0 harris-index=on color-survey=on hue=Orientation sat=Coherency bri=Original-Image");
	
	selectWindow("Color-survey-1");
	setBatchMode("show"); // Safety option
	run("Duplicate...", "tobesaved");
	saveAs("tiff", output + File.separator + "Cell" + File.separator + "Img_" + serial+1 + "c - Orientation - " + file);
	close(); // Close saved tiff file
	
	selectWindow("Color-survey-1");
	run("HSB Stack");
	run("Duplicate...", "title=dup-bright duplicate range=3-3"); // Use Brightness-slice -> Original Binarized Image as Mask
	run("8-bit"); // ...redundant
	run("Divide...", "value=255"); // {0,1} Values (Multiplicative Mask)
	
	selectWindow("Color-survey-1");
	run("Duplicate...", "title=dup-coher duplicate range=2-2"); // Use Saturation-slice -> Coherency
	run("8-bit"); // ...redundant
	imageCalculator("Multiply create", "dup-coher", "dup-bright"); // Mask Coherency to make weights for histogram
	selectWindow("Result of dup-coher");
	rename("weight");
	
	selectWindow("dup-bright");
	close();
	selectWindow("dup-coher");
	close();
	
	selectWindow("Color-survey-1");
	run("Duplicate...", "title=dup-orient duplicate range=1-1"); // Use Hue-slice -> Orientation
	run("8-bit"); // ...redundant
	setMinAndMax(0, 361); // Remap 256 (8-bit) to 180 levels -> rebin histogram! 361 has been empirically tested (256:180=x:256 -> x=364)
	run("Apply LUT");
	run("Add Specified Noise...", "standard=1"); // Noise will smooth histogram to avoid quantization spikes
	
	selectWindow("Color-survey-1");
	close();
	
	selectWindow("weight");
	getStatistics(areaW, meanW, minW, maxW);
	step = floor(maxW/10); // Ten-level weights (10 instead of 256 to speed up the computation)
	//print("Max of Coherency = " + maxW); // Just for testing purpose
	//print("Ten-level step = " + step); // Just for testing purpose
	
	mainhisto = newArray(180);

	if(step >= 1) { // ...that is (maxW < 10)
		for(h=0; h < 10; h++) {
			
			// Just for testing purpose
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
			
			setOption("BlackBackground", false);
			run("Convert to Mask");
			run("Divide...", "value=255"); // {0,1} Values (Multiplicative Mask)
			
			imageCalculator("Multiply create", "temp", "dup-orient");
			selectWindow("Result of temp");
			getStatistics(areaNoUse, meanNoUse, minNoUse, maxNoUse, stdNoUse, histo);
			//getHistogram(0, counts, 256); // As an alternative, to return just the histogram
			
			selectWindow("Result of temp");
			close();
			selectWindow("temp");
			close();
			
			histo[0] = 0; // Delete background bin (and zero-Coherency pixels)
			for(k=0; k<180; k++){
				mainhisto[k] = mainhisto[k] + histo[k]*((h+1)*step-step/2); // Build weighted histogram
			}
		}
	}
	else{
		//print("unique level", ":", 1, "-", maxW); // Just for testing purpose
		
		selectWindow("weight");
		run("Duplicate...", "title=temp");		
		setAutoThreshold("Default dark");
		call("ij.plugin.frame.ThresholdAdjuster.setMode", "B&W");
		setThreshold(1, maxW);
		setOption("BlackBackground", false);
		run("Convert to Mask");
		run("Divide...", "value=255"); // {0,1} Values (Multiplicative Mask)
		imageCalculator("Multiply create", "temp", "dup-orient");
		selectWindow("Result of temp");
		getStatistics(areaNoUse, meanNoUse, minNoUse, maxNoUse, stdNoUse, histo);
		
		selectWindow("Result of temp");
		close();
		selectWindow("temp");
		close();
		
		histo[0] = 0; // Delete background bin (and zero-Coherency pixels)
		for(k=0; k<180; k++){
			mainhisto[k] = histo[k]*4.5; // Build weighted histogram (4.5 is the midpoint of the coherency range [1-9])
		}
	}
	
	selectWindow("dup-orient");
	close();
	selectWindow("weight");
	close();
	
	mainhisto[0] = round((mainhisto[1] + mainhisto[179])/2); // Circular junction
	
	if(full) {
		Plot.create("MainHistogram", "Value", "Count", mainhisto);
		Plot.show();
		selectWindow("MainHistogram");
		saveAs("tiff", output + File.separator + "Cell" + File.separator + "Histo_" + serial+1);
		close();
	}
	
	selectWindow("MegaMatrix2");
	for (row = 0; row < 180; row++) {
		setPixel(row, serial+2, mainhisto[row]);
	}
	
	//return something;
}

// Plot MegaMatrix2
function PlotMM2(sampleNum) {
	
	selectWindow("MegaMatrix2");
	// LUTs can't be applied to 32-bit images -> Convert to 16-bit
	run("Conversions...", "scale"); // Scale from min–max to 0–65535 when converting from 32–bit to 16–bit
	getStatistics(areaM, meanM, minM, maxM);
	setMinAndMax(0, maxM);
	run("16-bit");
	
	makeRectangle(2, 0, sampleNum, 180); // Enhance Contrast to normalize Soma and Cytoskeleton separately, allowing for comparison
	run("Enhance Contrast", "saturated=0.35");
	run("Apply LUT");
	makeRectangle(sampleNum+3, 0, sampleNum, 180);
	run("Enhance Contrast", "saturated=0.35");
	run("Apply LUT");
	
	makeRectangle(2, 0, 2*sampleNum+1, 180); // Crop "Degree" column
	run("Crop");
	
	// Draw a separation white line between Shape and Cytoskeleton Orientations and apply a LUT
	getStatistics(areaMM, meanMM, minMM, maxMM);
	for (row = 0; row < 180; row++) {
		setPixel(sampleNum, row, maxMM);
	}
	run("Fire");
	//run("Apply LUT"); // Choose whether to apply it or not
	
	newHeight = 4*getHeight(); // Magnification 4x
	newWidth = 4*getWidth();
	run("Size...", "width=" + newWidth + " height=" + newHeight + " constrain average interpolation=None");
	run("Canvas Size...", "width=" + getWidth()+60 + " position=Center-Right zero");
	
	// Draw degree vertical axis
	newImage("VerticalAxis", "16-bit black", 60, 720, 1);
	setColor("white");
	setFont("SansSerif", 16, "bold");
	setJustification("right");
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
	
	makeRectangle(0, 0, 60, 720); // Crop "Degree" column
	run("Copy");
	selectWindow("MegaMatrix2");
	makeRectangle(0, 0, 60, 720);
	run("Paste");
	selectWindow("VerticalAxis");
	close();
	selectWindow("MegaMatrix2");
	run("Select None");

	//return something;
}
