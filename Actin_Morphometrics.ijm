macro "Actin Morphometrics" {

///////////////////////////////////////////////////////////////////////////////

// PURPOSE
// detect and analyze actin filaments
// processes all TIFF files in a directory

// INPUTS
// - asks user for a directory (folder) with all images that should be analyzed
// - images have to be in 8-bit TIFF format
// - all images must have their resolution (pixel size in micrometers) defined
// -- alternatively, a common resolution can be defined at runtime

// OUTPUTS
// for each input image:
//   a smoothed and background subtracted version of each input image ("sb")
//   a convoluted image with filaments highlighted ("sbf")
//   a "Mask" image that traces detected actin filaments for each input image
//   a "Filaments" image that shows actin filament intensity along the mask (with background = 0)
//   a "Distance Map" (32-bit) image that shows distance to nearest actin filament
//   an "Angles" image that color codes the angle of filaments (relative to horizontal)
//   a results file with various parameters (explanation in separate file)
// Summary table with parameters for all images in the directory
//   (user selects whether to report just a few or all parameters)
//   (a second copy of the table also includes simple statistics of the parameters)

// PROCEDURE
// - user sets a number of parameters
// - user selects directory (folder) that contains images
// -- macro detects actin filaments in each image and saves results to files


// VERSIONS
// 1  - AN.150130 - first implementation
// 2  - AN.150201 - also save info about individual filaments
// 3  - AN.150317 - include angle calculation (in internal function)
// 7  - AN.201019
// 8  - AN.201231 - streamline convolution and calculations
// 9  - AN.210115 - complete rewrite with many different function and new calculations
// 10 - AN.210217 - a few small edits to fix typos and tidy things up a bit
//                - corrected skewness calculation to include saturated pixels
//                - include calculation of black and saturated pixels
// 10a- AN.210123 - also determine detection threshold
//                - offer option to save all or just a few columns
// 10b- AN.210123 - also determine 75th and 90th percentile of filament intensities
// 10c- AN.210310 - also calculate Branch Density (branches per micron)
// 11 - AN.210318 - changed thresholding from upper 10 % to mean of pixels > 1 (auto)
//                - cleaned up and corrected several calculations
// 11a- AN.210319 - percentile calculation for filaments instead of cell
// 11b- AN.210405 - allow input of a manual ROI for cell area
// 12b- AN.230210 - spatial calibration can be entered manually for images
//                - analysis of "holes" without actin filaments (fixed threshold)
// 12c- AN.230410 - allow user to define distance between holes and filamnts
// 12d- AN.230808 - now can calculate holes info for entire image (or cell area only)
// 13 - AN.230905 - add absolute angle calculation relative to horizontal (+ variation)
// 13a- AN.230906 - also calculate angle variation (sd) and order parameter (align)
// 13b- AN.230907 - include output of weighted angles etc
// 14 - AN.230911 - move angle calculations to separate function + add angle image
// 14a- AN.230914 - added array to hold filament angle distribution
// 14b- AN.230919 - added order parameter for filament orientation + better angle calculation
// 14c- AN.230920 - also calculate bundling parameter (CV/FD)
// 15 - AN.230921 - offer option to do calculation for full frame of for convex hull (= cell area)
//                - save filament intensity histogram in Results file per image
//                - clean-up: report only a few parameters, remove strange ones (px intensities, area analysis,...)
//                - remove unnecessary calculations (holes...)
// 15a- AN.230929 - added quality check for too many saturated pixels (1% cutoff)
// 15b- AN.240105 - replaced thresholding method: mean (pixels > 0) -> Otsu (all pixels)
// 15c- AN.240110 - removed Otsu threshold and added Bundling=CV/Occ calculation
// 16 - AN.240419 - calculate Order Parameter relative to Mean Angle.
// 16a- AN.240424 - calculate Order Parameter relative to Mean Weighted Angle.
//                - calculate angles also as absolute deviation from horizontal.
// 16a1-AN.240425 - reversed up/down test of angles
// 16a2-AN.240426 - alternating angle identification
// 17 - AN.240429 - revised angle calculation to pick continuation with smallest turn    
// 17a- AN.240507 - fixed weighted Order Parameter calculation 
// 17b- AN.240508 - report Order Parameter relative to horizontal and mean angle
// 18 - AN.240607 - calculate Bundle Parameter as CV x median distance 
//                - changed order of functions to accommodate previous calculation
//                - removed tables and replaced with arrays
//                - reduced parameter list is now: 
//                    Occupancy, Med.Distance, Bundle Parameter, Mean Weighted Angle, Angular Variation, Order Parameter
// 18p- AN.240725 - cleanup of code to streamline calculations and output for publication.
///////////////////////////////////////////////////////////////////////////////

// // // general checks and settings
if (getVersion() < "1.54g") setOption("ExpandableArrays",true);
run("Set Measurements...", "decimal=5");  // report all measurements with 5 decimal places

// definitions:
tab = "\t";  
pixelIntensities = newArray(256);         // pixel intensities along mid-line of filaments
medDist = 0;                              // initial value that will be changed later


param_info = newArray(30); // to hold the names of the measured parameters (first element is empty)
param = newArray(30);      // to hold the values of the measured parameters (first element is number of parameters)
param_main = newArray(10); // to hold the positions of those parameters that will be in final table
                           //   (first element is number of parameters)
// these arrays will be passed into every function >> param_info, param, param_reported 

// default values:
smoothingRadius = 0.2;              // radius for image smoothing (in microns)
kernelSize = 2.5 * smoothingRadius; // length of convolution kernels (in microns)
pixelSize = 0;                      // pixel size (only used for images that are not pre-calibrated)
smoothingAngles = 0.5;              // length of lines for calculating angles (in microns)

fullFrame = 0;  // set this to 1, if calculations should be based on the entire image frame
outputAll = 0;  // set this to 1, if only all parameters should be reported in the summary table



// // // ask user for confirmation of default values

Dialog.create("Parameter Definition");
Dialog.addNumber("Smoothing radius (in microns) ",smoothingRadius);
Dialog.addNumber("Linear enhancement kernel (in microns)",kernelSize);
Dialog.addNumber("Filament smoothing for angle measurements (in microns)",smoothingAngles);
Dialog.addNumber("Pixel size for uncalibrated images (in microns)", pixelSize);
Dialog.addCheckbox("Report density based on full image frame",fullFrame);
Dialog.addCheckbox("Report all results in summary table",outputAll);
Dialog.show();

smoothingRadius = Dialog.getNumber();
kernelSize = Dialog.getNumber();
smoothingAngles = Dialog.getNumber();
pixelSize = Dialog.getNumber();
fullFrame = Dialog.getCheckbox();
outputAll = Dialog.getCheckbox();


// // // ask user for directory with images

showStatus("Select the directory that contains the images:");
path = getDirectory("Choose a Directory");
list = getFileList(path);

setBatchMode(1);

// create table to collect results from individual images into a large table
allTable = "Measured Parameters";
Table.create(allTable);            // this table will contain the reported results
allTablePlusStats = "Measured Parameters with Statistics";
Table.create(allTablePlusStats);   // same as before, but with mean, sd, min, max, etc


// process all files //////////////////////////////////////////////////

nI = 0;                             // keep track of number of images
for (f=0; f<list.length; f++) {     // for each file in directory ... ***
  if (endsWith(list[f],".tif")) {   // only work with TIFF files... ==
    nI++;
    open(list[f]);                  // open image file
    image = getImageID();           // get image ID of original image
    run("8-bit");                   // make sure this is processed as an 8-bit image
    showStatus("processing file "+list[f]);
    extension = lastIndexOf(list[f],".tif"); // find where extension starts
    imageName = substring(list[f],0,extension); // remove extension 
    getPixelSize(unit,pxW,pxH);
    if (pxW != 1) { // pixel dimension defined?
      pxSize = (pxW+pxH)/2; // this is needed since the confocal produces not quite square pixels?
    } else { // pixel dimension not defined
      pxSize = pixelSize;
      setVoxelSize(pxSize,pxSize,1,"micron");
    } // END IF pixel dimensions defined?
    // define/reset arrays for measured parameters
    param_info[0] = "Parameter Names"; // info
    param[0] = 0;                      // so far no parameters measured
    param_main[0] = 0;                 // so far no main parameters measured

    roiManager("reset");              // clear existing list of selections
    if (selectionType >= 0) { // is the cell area already selected?
      setSelectionName("CellArea");   // is this needed <<<<<<<<<<
      roiManager("add");              // add this selection to ROI manager
      cellROI = roiManager("count")-1;  // remember the ROI
      saveAs("selection",path+imageName+"-cell.roi"); // also save as separate file
    } // cell area already selected?

    // call processing and analysis functions: ////////////////////////////

    smBkgSubImage = pre_process(image, smoothingRadius);
        showStatus("processing file "+list[f]+"  linear enhancements");
    filamentsImage = linear_convolutions(smBkgSubImage, kernelSize);
       showStatus("processing file "+list[f]+"  skeletonization");
    skeletonImage = skeleton(filamentsImage);
        showStatus("processing file "+list[f]+"  measuring density");
    filament_properties(skeletonImage, fullFrame, param_info, param, param_main);
        showStatus("processing file "+list[f]+"  measuring distances");
    distMapImage = distances(skeletonImage, fullFrame, param_info, param, param_main);
    for (i=1; i<param[0]+1; i++) { // FOR find the Median Distance value
      if (param_info[i] == "Median Distance") { medDist = param[i]; }
    } // END FOR find the Median Distance value
        showStatus("processing file "+list[f]+"  measuring bundling levels");
    intensityImage = bundling(image, skeletonImage, medDist, param_info, param, param_main, pixelIntensities);
        showStatus("processing file "+list[f]+"  measuring angles and ordering");
    angleImage = angles(intensityImage, smoothingAngles, param_info, param, param_main);
        showStatus("processing file "+list[f]+"  final processing");


    // save results: ////////////////////////////////////////////////////////////

    // 1. Save all information for this image in a separate file.
    out = File.open(path+imageName+"-Results.txt"); // open text file for all results
    print(out,"Image Name:\t"+list[f]);
    print(out,"Pixel Size\t"+pxSize);
    print(out,"");
    for (i=1; i<param[0]+1; i++) { // FOR all parameters
      print(out,param_info[i]+tab+param[i]);
    } //END FOR all parameters
    print(out,"");
    print(out,"pixel intensity distribution along filaments");
    print(out,"");
    print(out,"intensity"+tab+"count");
    for (i=0; i<256; i++) { print(out,i+tab+pixelIntensities[i]); }
    print(out,"");


    // 2. Save (all or select) parameters into table
    Table.set("File",nI-1, list[f], allTable);         // write file name into results table
    if (outputAll) { // IF save all parameters ###
      Table.set("Pixel Size",nI-1, pxSize, allTable);  // write pixel size into summary table
      for (i=1; i<param[0]+1; i++) { // FOR all parameters *
        Table.set(param_info[i], nI-1, param[i], allTable);
      } //END FOR all parameters *
    } else { // ELSE save select parameters ###
      for (i=1; i<param_main[0]+1; i++) { // FOR all parameters in the param_main list **
        p = param_main[i];              // get the number of the parameter
        Table.set(param_info[p], nI-1, param[p], allTable);
      } // END FOR all parameters in the param_main list **
    } // END IF save all or select parameters ###

    // clean up after processing the file ////////////////////////////////////////////
    File.close(out); // close results file

    // save images and close them
    selectImage(smBkgSubImage);
    name = getTitle();
    saveAs("tiff",path+name);
    selectImage(filamentsImage);
    name = getTitle();
    saveAs("tiff",path+name);
    selectImage(skeletonImage);
    name = getTitle();
    saveAs("tiff",path+name);
    selectImage(intensityImage);
    name = getTitle();
    saveAs("tiff",path+name);
    selName = imageName+"-Filaments.roi";
    saveAs("selection",path+selName);
    selectImage(angleImage);
    name = getTitle();
    saveAs("tiff",path+name);
    selectImage(distMapImage);
    name = getTitle();
    saveAs("tiff",path+name);
    selName = imageName+"-CellArea.roi";
    saveAs("selection",path+selName);
    close("*"); // close all image windows
  } // END IF .. tiff file ==
} // END FOR .. all files ***


// general cleanup at the very end \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

imagesN = table_stats(allTable, allTablePlusStats);
Table.save(path+allTable+".txt",allTable);
Table.save(path+allTablePlusStats+".txt",allTablePlusStats);
close(allTable);


setBatchMode(0);
updateDisplay();
beep();
showStatus("Done! Processed "+toString(imagesN)+" images.");

///////////////////////////////////////////////////////////////////////////////

} // END MACRO "actin morphometrics"

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////

function pre_process(imageID, radius) {

///////////////////////////////////////////////////////////////////////////////

// PURPOSE:
// prepare image for detection of actin filaments

// INPUT
// imageID = pointer to image window
// radius = radius (in microns) for smoothing

// RETURN VALUE
// pointer to modified image

// PROCEDURE
// 1. get image parameters: pixel size, size, name
// 2. create a copy of the image for processing
// 3. calculate smoothing radius and rolling ball radius in pixel
//    (rolling ball is 10-fold larger)
// 4. perform smoothing and background subtraction

// VERSIONS
// 1.0 = AN.201230

///////////////////////////////////////////////////////////////////////////////

// 1. get image parameters: pixel size, title, size
selectImage(imageID);
getPixelSize(unit,pxW,pxH);
px = pxW;
w = getWidth();
h = getHeight();
original = getTitle();
extension = lastIndexOf(original,".tif"); // find where extension starts
basename = substring(original,0,extension); // remove extension 
newName = basename + "-sb.tif";   // name of the processed image will end in "sb"

// 2. create a copy of the image for processing
run("Select All"); run("Copy");   // copy image data
newImage(newName,"8-bit",w,h,1);  // 
setVoxelSize(px,px,1,"micron");   // define pixel size
run("Paste");                     // paste image data into new window
processed = getImageID();

// 3. calculate smoothing radius and rolling ball radius in pixel
rSm=radius/px; // radius for smoothing in pixels
rBkg = 10*rSm; // radius for rolling ball background subtraction

// 4. perform smoothing and background subtraction
run("Mean...", "radius="+rSm);
run("Subtract Background...", "rolling="+rBkg+" sliding disable");

// 5. cleanup
return processed;

///////////////////////////////////////////////////////////////////////////////

} // END of FUNCTION pre_process

///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////

function linear_convolutions(imageID, sizeK) {

///////////////////////////////////////////////////////////////////////////////

// PURPOSE:
// enhance linear structures such as actin filaments

// INPUT
// imageID = pointer to image window
// sizeK = kernel size (in microns)

// RETURN VALUE
// pointer to modified image

// PROCEDURE
// 1. get image parameters: pixel size, size, name
// 2. create four copies of the image in a stack for processing
// 3. calculate kernel size in pixel
// 4. run convolution with different kernel on each slice of stack
//    four directions: horizontal, vertical, down diagonal, up diagonal
// 5. combine four slices into one
// 6. remove stray spots outside of cell (if pre-selected)

// VERSIONS
// 1.0 = AN.201231
// 2.0 = AN.210405 
//     - added removal of signal outside of (manually) selected cell area

///////////////////////////////////////////////////////////////////////////////

// 1. get image parameters: pixel size, title, size
selectImage(imageID);
getPixelSize(unit,pxW,pxH);
px = (pxW+pxH)/2; // this is needed since the confocal produces not quite square pixels?
w = getWidth();
h = getHeight();
original = getTitle();
extension = lastIndexOf(original,".tif");   // find where extension starts
basename = substring(original,0,extension); // remove extension 

// 2. create four copies of the image in a stack for processing
run("Select All"); run("Copy");             // copy image data
newImage("tempStack","8-bit",w,h,4);        // create temporary stack for convolutions
for (i=1; i<=4; i++) { 
  setSlice(i);                              // for each slice
  run("Paste");                             // fill with starting image
  } //END FOR  each slice
stack = getImageID();

// 3. calculate kernel size in pixel
kernelSize = round(sizeK/px);               // first approximation of kernel size
if (kernelSize%2 == 0) { kernelSize++; }    // if even, increment size
if (kernelSize < 5) { kernelSize = 5; }     // if too small, increase to 5 pixels

// 4. run convolution with different kernel on each slice of stack
setSlice(1);
parameter = "text1=["+kernel_line_horiz(kernelSize)+"] normalize slice";
run("Convolve...", parameter);
setSlice(2);
parameter = "text1=["+kernel_line_vert(kernelSize)+"] normalize slice";
run("Convolve...", parameter);
setSlice(3);
parameter = "text1=["+kernel_line_down(kernelSize)+"] normalize slice";
run("Convolve...", parameter);
setSlice(4);
parameter = "text1=["+kernel_line_up(kernelSize)+"] normalize slice";
run("Convolve...", parameter);

// 5. combine four slices into one
run("Z Project...", "projection=[Max Intensity]");   // flatten stack
rename(basename+"f.tif");                            // add "f" to filename (for filaments)
setVoxelSize(px,px,1,"micron");                      // define pixel size
convolved = getImageID();

// 6. remove everything outside the cell selection
if (RoiManager.size>0) {         // if cell area was already selected
  roiManager("select",cellROI);  // restore cell selection
  run("Make Inverse");           // select outside of cell
  setColor(0);                   // paint in black
  fill();
} // cell are already selected
  

// 6. cleanup
close("temp*");   // close temporary stack
return convolved; // return imageID of convolved image

///////////////////////////////////////////////////////////////////////////////

} // END of FUNCTION linear_convolutions

///////////////////////////////////////////////////////////////////////////////




// ** functions that define the kernel strings for linear convolution **

////////////////////////////////////////////////////////////////////////////////////////////////

function kernel_line_horiz(size) {

// PURPOSE:
// defines a string that can be used for convolution command
// kernel enhances linear features that run from top to bottom
// flanking pixels are weighted -1, central pixel is sum of flanking
// parameter defines the length of the string.

// RETURN VALUE:
// string that defines the convolution kernel


if (size%2 == 0) { size++; } // make sure the size of the kernel is an odd number

kernel = ""; // start with an empty kernel
for (i=0; i<((size-1)/2); i++) { // for the first half of positions
  kernel+="-1 ";
  } // END FOR first half
kernel+=toString(size-1)+" "; 
for (i=0; i<((size-1)/2); i++) { // for the second half of positions
  kernel+="-1 ";
  } // END FOR second half
kernel+="\n";

return kernel;

} // END FUNCTION

////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////

function kernel_line_vert(size) {

// PURPOSE:
// defines a string that can be used for convolution command
// kernel enhances linear features that run from left to to right
// flanking pixels are weighted -1, central pixel is sum of flanking
// parameter defines the length of the string.

// RETURN VALUE:
// string that defines the convolution kernel


if (size%2 == 0) { size++; } // make sure the size of the kernel is an odd number

kernel = ""; // start with an empty kernel
for (i=0; i<((size-1)/2); i++) { // for the first half of positions
  kernel+="-1\n";
  } // END FOR first half
kernel+=toString(size-1)+"\n"; 
for (i=0; i<((size-1)/2); i++) { // for the second half of positions
  kernel+="-1\n";
  } // END FOR second half

return kernel;

} // END FUNCTION

////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////

function kernel_line_down(size) {

// PURPOSE:
// defines a string that can be used for convolution command
// kernel enhances linear features that run from lower left to upper right
// flanking pixels are weighted -1, central pixel is sum of flanking
// parameter defines the length of the string.

// RETURN VALUE:
// string that defines the convolution kernel

if (size%2 == 0) { size++; } // make sure the size of the kernel is an odd number

kernel = ""; // start with an empty kernel
for (i=0; i<size; i++) { // go through all lines
  for (j=0; j<size; j++) { // fill all positions in line
    if (j == i) { // on diagonal?
      if ((i == (size-1)/2)) { // in center line?
        kernel+=toString(size-1)+" "; 
        } else {
        kernel+="-1 "; 
        } // END IF center spot
      } else { 
      kernel+="0 "; 
      } // END IF on diagonal
    } // END FOR all positions
  kernel+="\n";
  } // END FOR all lines

return kernel;

} // END FUNCTION

////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////

function kernel_line_up(size) {

// PURPOSE:
// defines a string that can be used for convolution command
// kernel enhances linear features that run from upper left to lower right
// flanking pixels are weighted -1, central pixel is sum of flanking
// parameter defines the length of the string.

// RETURN VALUE:
// string that defines the convolution kernel

if (size%2 == 0) { size++; } // make sure the size of the kernel is an odd number

kernel = ""; // start with an empty kernel
for (i=0; i<size; i++) { // go through all lines
  for (j=0; j<size; j++) { // fill all positions in line
    if (j == (size-i-1)) { // on diagonal?
      if ((i == (size-1)/2)) { // in center line?
        kernel+=toString(size-1)+" "; 
        } else {
        kernel+="-1 "; 
        } // END IF center spot
      } else { 
      kernel+="0 "; 
      } // END IF on diagonal
    } // END FOR all positions
  kernel+="\n";
  } // END FOR all lines

return kernel;

} // END FUNCTION

////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////

function skeleton(imageID) {

///////////////////////////////////////////////////////////////////////////////

// PURPOSE:
// enhance linear structures such as actin filaments

// INPUT
// imageID = pointer to image window
// fraction = fraction of pixels that should be considered signal

// RETURN VALUE
// pointer to modified image

// PROCEDURE
// 1. get image parameters: pixel size, size, name
// 2. create a copy of the image for processing
// 3. make binary and remove strays
// 4. measure area of filaments & cell and calculate occupancy
// 5. skeletonize 

// VERSIONS
// 1.0 = AN.201231
// 1.1 = AN.210198
//     - added calculation of occupancy
// 2.0 = AN.210318
//     - removed calculations since they were not correct
// 2.1 = AN.230206
//     - small fixes

///////////////////////////////////////////////////////////////////////////////

// 1. get image parameters: pixel size, size, name
selectImage(imageID);
getPixelSize(unit,pxW,pxH);
px = (pxW+pxH)/2; // this is needed since the confocal produces not quite square pixels?
w = getWidth();
h = getHeight();
original = getTitle();
extension = lastIndexOf(original,".tif");   // find where extension starts
basename = substring(original,0,extension); // remove extension 
newName = basename + "-Mask.tif";           // name of the processed image will end in "Mask"

// 2. create a copy of the image for processing
run("Select All"); run("Copy");   // copy image data
newImage(newName,"8-bit",w,h,1);  // 
setVoxelSize(px,px,1,"micron");   // define pixel size
run("Paste");                     // paste image data into new window
processed = getImageID();

// 3. make binary and remove strays
// setThreshold(t, 255);//
run("Options...", "iterations=1 count=1 black"); // set background to black
setForegroundColor(255,255,255);                 // thresholded area should be white
setAutoThreshold("Percentile dark no-reset");    // first identify all non-zero pixels
run("Create Selection");                         // use only non-zero pixels
setAutoThreshold("Mean dark no-reset");          // define threshold 
run("Convert to Mask");                          // turn into binary image
run("Open");                                 // to fill in small holes
run("Close-");                               // to remove single pixels

// 4. skeletonize 
run("Skeletonize");

// 5. cleanup
return processed; // return imageID of skeletonized image

///////////////////////////////////////////////////////////////////////////////

} // END of FUNCTION skeleton

///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////

function filament_properties(imageID, fullFrame, p_names, p_values, p_core) {

///////////////////////////////////////////////////////////////////////////////

// PURPOSE:
// calculates bundling parameters for image based on skeletonized mask

// INPUT
// imageID = pointer to image window with intensity modulated filaments
// fullFrame = flag for full frame (1) or cell area (0)
// p_names = pointer to array that holds parameter names
// p_values = point to array that holds parameter values
// p_core = point to array that identifies the position of key parameters

// RETURN VALUES
// pixels, total length (microns), filament density (microns per cell area), occupancy, 
// nBranches, branch freq, mean angle, parallelness
//   all returned in p_values array
//   core parameters are: occupancy
//     note: mean angle and parallelness calculated according to 
//             Ueda et al., (2010) PNAS 017: 6894-6899

// PROCEDURE
// 1. get image parameters: pixel size, size, name
// 2. initialize arrays etc
// 3. calculate average angles, paralleness, and branch points
// 4. determine occupancy and filament density
// 5. measure average angle and SD.
// *. write everything into parameterArrays

// VERSIONS
// 1.0 = AN.201231
// 1.1 = AN.210221
//     - corrected filament length output to microns
// 1.2 = AN.210310
//     - added Branch Density (number of branches per micron filament)
// 1.3 = AN.210318
//     - added calculation of occupancy and filament density
// 1.4 = AN.210405
//     - added possibility to have predefined cell area
// 2.0 = AN240613
//     - removed results table and replaced with arrays

///////////////////////////////////////////////////////////////////////////////

// 1. get image parameters: pixel size, title, size
selectImage(imageID);
getPixelSize(unit,pxW,pxH);
px = (pxW+pxH)/2; // this is needed since the confocal produces not quite square pixels?
w = getWidth();
h = getHeight();
original = getTitle();
extension = lastIndexOf(original,".tif"); // find where extension starts
basename = substring(original,0,extension); // remove extension 

// 2. initialize arrays etc
nBranches = 0; // no branch points so far
nPairs = newArray(4); // will hold the number of pixel pairs in four directions
    // 0 = right = 0°; 1 = down&right = 45°; 2 = down = 90°; 3 = down&left = 135°
for (i=0; i<4; i++) { nPairs[i]=0; } // no pixel pairs so far

// 3. calculate average angles, paralleness, and branch points
for (x=2; x<w; x++) {
  for (y=2; y<h; y++) {
    if (getPixel(x,y)>0) { // if pixel value > 0...
      hits=0;
      if (getPixel(x+1,y)>0) { // pixel to right continues
        nPairs[0]++;
        hits++;
      } // to right
      if (getPixel(x+1,y+1)>0) { // pixel to lower right continues
        nPairs[1]++;
        hits++;
      } // to lower right
      if (getPixel(x,y+1)>0) { // pixel to bottom continues
        nPairs[2]++;
        hits++;
      } // to right
      if (getPixel(x-1,y+1)>0) { // pixel to below left continues
        nPairs[3]++;
        hits++;
      } // to below left
      if (hits>1) {
        nBranches++;
      } // one more branch point
    } // END IF ... pixel value > 0 
  } // ENDFOR all y's
} // ENDFOR all x's

total = nPairs[0] + nPairs[1] + nPairs[2] + nPairs[3];   // all pixel pairs
length = nPairs[0] + nPairs[2] + sqrt(2)*(nPairs[1] + nPairs[3]); // add lengths of all pairs
length = length * px; // convert total filament length from pixels to microns
branchFreq = nBranches/total;                            // % branches
branchDens = nBranches/length;                           // branch density (per micron)

// calculate mean angle and parallelness according to Ueda 2010 PNAS
deltaDiag = abs(nPairs[1]-nPairs[3]);
deltaVert = abs(nPairs[0]-nPairs[2]);
at = atan(deltaDiag/(deltaVert+deltaDiag));
if (nPairs[0]>=nPairs[2]) {
  if (nPairs[1]>=nPairs[3]) {
    meanAngle = at;
  } else {
    meanAngle = PI - at;
  }
} else {
  if (nPairs[1]>=nPairs[3]) {
    meanAngle = PI/2 + at;
  } else {
    meanAngle = PI/2 - at;
  }
}
meanAngle = meanAngle * 180 / PI; // covert to degrees
if (meanAngle>90) { meanAngle -= 180; } // bring angles to between -90 and +90 degrees
parallelness = (deltaVert + deltaDiag) / total;


// 4. measure area of filaments & cell and calculate occupancy
run("Set Measurements...", "area redirect=None decimal=3"); // only measure area
setThreshold(1,255); 
run("Create Selection"); // turn threshold into selection of filaments
run("Measure");          // measure filament area

if (RoiManager.size>0) {         // if cell area was already selected
  roiManager("select",cellROI);  // restore cell selection
} else { // ... cell are not already selected
  run("Convex Hull");      // select area around all filaments
} // select cell area
run("Measure");          // measure cell area
run("Select None");      // clear selection
resetThreshold();         // turn of threshold

areaFilaments = getResult("Area",0);
areaCell = getResult("Area",1);
occupancyCA = areaFilaments/areaCell;      // occupancy relative to cell area 
occupancyFF = areaFilaments/(w*h*px*px);   // occupancy relative to full image frame
run("Clear Results");
filamentDensityCA = length/areaCell;       // filament density in cell area (micron / micron^2)
filamentDensityFF = length/(w*h*px*px);    // filament denisty in full image frame

if (fullFrame) { // IF full frame
  occupancy = occupancyFF;
} else { // ELSE cell area
  occupancy = occupancyCA;
} // END IF full frame or cell area

// 4. write everything into parameter arrays
p_count = p_values[0]; // number of parameters measured so far

p_names[p_count+1] = "Total Pixel Pairs";
p_names[p_count+2] = "Total Filament Length";
p_names[p_count+3] = "Cell Area";
p_names[p_count+4] = "Filament Density";
p_names[p_count+5] = "Occupancy";
p_names[p_count+6] = "# Branches";
p_names[p_count+7] = "Branch Freq.";
p_names[p_count+8] = "Branch Density";
p_names[p_count+9] = "Parallelness (Ueda(2010)";
p_names[p_count+10]= "Mean Angle (Ueda2010)";

p_values[p_count+1] = total;
p_values[p_count+2] = length;
p_values[p_count+3] = areaCell;
if (fullFrame) { // IF measure Density based on full image frame **
  p_values[p_count+4] = filamentDensityFF;
  p_values[p_count+5] = occupancyFF;
} else { // ELSE measure Density based on cell area **
  p_values[p_count+4] = filamentDensityCA;
  p_values[p_count+5] = occupancyCA;
} // END IF full frame or not? **
p_values[p_count+6] = nBranches;
p_values[p_count+7] = branchFreq;
p_values[p_count+8] = branchDens;
p_values[p_count+9] = parallelness;
p_values[p_count+10]= meanAngle;

p_values[0] = p_count+10; // 10 new parameters measured

p_core[0]++;
p_core[p_core[0]] = 5; // Occupancy is a core measurement

return;

///////////////////////////////////////////////////////////////////////////////

} // END of FUNCTION filament_properties

///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////

function bundling(imageIDoriginal, imageIDmask, medDist, p_names, p_values, p_core, intensityArray) {

///////////////////////////////////////////////////////////////////////////////

// PURPOSE:
// calculates bundling parameters for image based on skeletonized mask

// INPUT
// imageIDoriginal = pointer to original, unprocessed image window
// imageIDmask = pointer to image window with skeletonized filaments
// medDist = median distance to nearest filament
// p_names = pointer to array that holds parameter names
// p_values = point to array that holds parameter values
// p_core = point to array that identifies the position of key parameters
// parameterArray = pointer to array for distribution of pixel intensities

// RETURN VALUE
// pointer to intensity modulated filament skeleton
// mean intensity, sd, skewness, coefficient of variation, bundle parameter 
//    (these are passed back in p_values array)
//    (core parameters: bundle parameter (CV*distance)
// distribution of pixel intensities in filaments (stored in intensityArray)

// PROCEDURE
// 1. get image parameters: pixel size, size, name
// 2. use mask image to extract signal intensities along filaments
// 3. calculate skewness and coefficient of variation of pixel intensities

// VERSIONS
// 1.0 = AN.210105
// 1.1 = AN.210108
//     - removed occupancy calculation
// 1.2 = AN.210222
//     - included calculation of black and saturated pixels based on original image 
//     - determine dimmest pixel = threshold of detection
// 1.3 = AN.210226
//     - include calculation of 75th and 90th percentile of pixel intensities
// 1.4 = AN.210318
//     - removed some unnecessary calculations
// 1.5 = AN.210319
//     - changed 75percentile calculation to based on filaments instead of cell
// 1.6 = AN.210405
//     - added possibility to use predefined cell area
// 1.7 = AN.210920
//     - added calculation of bundling = coefficient of variation / filament density
// 2.0 = AN.210921
//     - introduced array with pixel intensities
//     - cleaning up code, removing some calculation, ...
// 2.1 = AN.240610
//     - calculating bundling as coefficient of variation * median distance
// 3.0 = AN.240613
//     - removed results table and replaced with results array

///////////////////////////////////////////////////////////////////////////////

// 1. get image parameters: pixel size, size, name, cell area

selectImage(imageIDmask);
getPixelSize(unit,pxW,pxH);
px = (pxW+pxH)/2; // this is needed since the confocal produces not quite square pixels?
maskName = getTitle();
extension = lastIndexOf(maskName,"-Mask.tif"); // find where extension starts
basename = substring(maskName,0,extension);    // remove extension 

// 2. use mask image to extract signal intensities along filaments
imageCalculator("and create",imageIDoriginal,imageIDmask);
rename(basename+"-filaments.tif"); // intensity modulated filaments
intensityModulated = getImageID;
setVoxelSize(px,px,1,"micron");   // define pixel size


// 3. calculate skewness and coefficient of variation of pixel intensities

setThreshold(1,255);     // find filaments and ignore background
run("Create Selection"); // select all filaments
run("Set Measurements...", "area mean standard skewness redirect=None decimal=3");
run("Measure");
filMean = getValue("Mean");
filSD = getValue("StdDev");
filSkew = getValue("Skew");
cov = filSD/filMean;
bundled = cov*medDist;
getHistogram(values,intensityArray,256);


// 4. save results

//// save values in parameter arrays

p_count = p_values[0]; // number of parameters measured so far

p_names[p_count+1] = "Mean Intensity [Filaments]";
p_names[p_count+2] = "StdDev Intensity [Filaments]";
p_names[p_count+3] = "Skewness";
p_names[p_count+4] = "Coefficient of Variation";
p_names[p_count+5] = "Bundle Parameter";

p_values[p_count+1] = filMean;
p_values[p_count+2] = filSD;
p_values[p_count+3] = filSkew;
p_values[p_count+4] = cov;
p_values[p_count+5] = bundled;

p_values[0] = p_count + 5;  // five more parameters measured

p_core[0]++; // one more core parameter measured
p_core[p_core[0]] = p_count+5; // Bundle Parameter is a core measurement



return intensityModulated; // return imageID of intensity modulated skeleton image


///////////////////////////////////////////////////////////////////////////////

} // END of FUNCTION bundling

///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////

function angles(intensityImage, smoothingLength, p_names, p_values, p_core) {

///////////////////////////////////////////////////////////////////////////////

// PURPOSE:
// calculates angles

// INPUT
// intensityImage = pointer to image window with intensity-modulate skeletonized filaments
// smoothingLength = length of line that will be used to measure angles of filaments (microns)
// p_names = pointer to array that holds parameter names
// p_values = point to array that holds parameter values
// p_core = point to array that identifies the position of key parameters

// RETURN VALUE
// pointer to image with filament angles color coded (angleImage)
// mean angle (*), angle standard deviation (*), mean wheighted angle, angular variation, 
//  order parameter (relative to horizontal), order parameter (realtive to mean weighted angle), 
//  weighted angle deviation
//     (*) according to Madison et al. (2015) Plant Physiol. 169: 1946-1960
//     (passed back in p_values)
//     (core parameters: mean weighted angle, angular variation, order parameter (relative to mean weighted angle)


// PROCEDURE
// 1. copy intensity image to new temporary window (+ create new empty anglesImage)
// 2. scan through all pixels to find a filament and copy coordinates into array
//    (also erase filament in temp window to prevent double counting)
// 3. measure angles along filament over larger distance (2*range)
// 4. calculate mean angle, sd angle, order parameter
// 5. plot angles in anglesImage
// 6. after all filaments are done, write all results to arrays

// VERSIONS
// 1.0 = AN.230912
//     - adapt from pollen actin macro (actin-longitudinal-v5a-AN.150704)
//     - added plot of angles in new image
// 1.1 = AN.230914
//     - added array to hold filament angles (in bins of 5° each)
// 2.0 = AN.230919
//     - added order parameter (how much are filaments aligned with horizontal?)
//          S = mean of [(2 * (cos(angle))^2) - 1]
//     - removed filament angle array
//     - direct calculation of mean/sd/order and weighted by brightness
// 3.0 = AN.240419
//     - calculating order parameter relative to mean angle (unweighted)
//     - calculating means/sds from 32-bit images 
// 3.1 = AN.240422
//     - calculating order parameter relative to mean angle (weighted)
//     - calculating mean deviation from horizontal as absolute value of weighted angles
// 3.2 = AN.240425
//     - picking continuation pixels alternates between up and down bias
// 4.0 = AN.240429
//     - picking continuation of filament with smallest turn angle to remove bias
// 4.1 = AN.240507
//     - fixed weighting calculation for order parameter
// 4.2 = AN.240508
//     - calculate order parameter relative to both horizontal and for mean angle
// 5.0 = AN.240613
//     - removed results table and replaced with results arrays

///////////////////////////////////////////////////////////////////////////////

// general parameters:
range = smoothingLength/2;   // distance of line ends from mid point for angle measurements

// 1. copy image to new window to be able to erase pixels.
selectImage(intensityImage);
w = getWidth(); h = getHeight();              // image dimensions (in pixels)
getPixelSize(unit,pxW,pxH); px = (pxW+pxH)/2; // pixel size (in "unit")
filamentsName = getTitle();
extension = lastIndexOf(filamentsName,"-filaments.tif"); // find where extension starts
basename = substring(filamentsName,0,extension);         // remove extension 
anglesName = basename+"-Angles.tif";                     // add new extension

run("Select All");
run("Copy");
newImage("temp", "8-bit black", w, h, 1);      // empty image to tick off filaments
run("Paste"); run("Select None");              // plot filaments in the temporary image
temporary = getImageID();                      // remmeber the handle
newImage(anglesName, "8-bit black", w, h, 1);  // to show filament angles
angleImage = getImageID();                     // remember the handle
newImage(anglesName, "32-bit black", w, h, 1);  // to hold filament angles as floating point numbers
run("Set...", "value=NaN");                    // set all angles to nonexistent
angleRadians = getImageID();                   // remember the handle
selectImage(temporary);    // all measurements will be done on the temporary image

// set up variables and arrays
filX = newArray(2*w);      // array to hold X coordinates of filament
filY = newArray(2*w);      // array to hold Y coordinates of filament
filBright = newArray(2*w); // array to hold brightness along filament
nFil = 0;                  // number of filaments
meanAng = 0; meanAngW = 0; // mean filament angle (horizontal = 0) (unweihghted and weighted)
ssqAng = 0; ssqAngW = 0;   // sum of squares of all angles (unweighted and weighted)
alignAng=0; alignAngW=0;   // order parameter for absolute angles (<2*cos^2(ang) - 1>)
nAng = 0;                  // number of angle measurements
tBright = 0;               // total brightness of filament
range = range / px;        // convert range from microns to pixels

// 2. scan image from left to right for filaments and write positions into array
for (xx=0; xx<w; xx++) {
  for (yy=0; yy<h; yy++) {
    if (getPixel(xx,yy)>0) { // if pixel value > 0 (i.e. if part of filament) ... continue analysis
      x = xx;  y = yy;                     // convert to local variables
      nFil++;                              // one more filament
      filLen = 0;                          // filament length (starts at 0)
      filX[filLen] = x; filY[filLen] = y;  // store pixel position
      filBright[filLen] = getPixel(x,y);   // store pixel brightness
      filLen++;                            // filament longer 
      lastAng = 0; preLastAng = 0;         // prefer horizontal filaments
      setPixel(x,y,0);                     // remove pixel to prevent double counting
      endFil = 0;       // filament not yet finished
      // find filament coordinates  // // // // // // // // // // // // // // //
      while (!endFil) { // continue until end of filament 
        endFil=1;       // assume filament will end with this pixel
        turn = 180;     // assume worst possible turn angle
        nextAng = 180;  // assume worst possible continuation
        priorAng = (lastAng + preLastAng)/2; // use average of the last two angles
        // check all 8 neighbors for next pixel in filament

        up = getPixel(x,y-1);     // get pixel intensity above current spot
        if (up>0) {               // if part of filament
          endFil = 0;             // filament continues
          newAng = -90; thisTurn = abs(priorAng) - abs(newAng);    // compare to old direction
          if (abs(thisTurn)<abs(turn)) {               // if less of a kink...
            // preLastAng=lastAng; lastAng=newAng; 
            bright=up; nextAng=newAng; turn=thisTurn; nextX=x;nextY=y-1; // save info
          } else if (abs(thisTurn)==abs(turn)) {       // if same kink as before ...
            if (abs(newAng)<abs(nextAng)) {            // only save if closer to horizontal
              // preLastAng=lastAng; lastAng=newAng; 
              bright=up; nextAng=newAng; turn=thisTurn; nextX=x;nextY=y-1; // save info
            } // END IF - more horizontal
          } // END IF - same kink as before
        } // END IF - part of filament

        upright = getPixel(x+1,y-1);   // get pixel intensity above&right of current spot
        if (upright>0) {               // if part of filament
          endFil = 0;   // filament continues
          newAng = -45; thisTurn = abs(priorAng) - abs(newAng);  // compare to old direction
          if (abs(thisTurn)<abs(turn)) {               // if less of a kink...
            // preLastAng=lastAng; lastAng=newAng; 
            bright=upright; nextAng=newAng; turn=thisTurn; nextX=x+1;nextY=y-1;// save info
          } else if (abs(thisTurn)==abs(turn)) {       // if same kink as before ...
            if (abs(newAng)<abs(nextAng)) {            // only save if closer to horizontal
              // preLastAng=lastAng; lastAng=newAng; 
              bright=upright; nextAng=newAng; turn=thisTurn; nextX=x+1;nextY=y-1; // save info
            } // END IF - more horizontal
          } // END IF - same kink as before
        } // END IF - part of filament

        right = getPixel(x+1,y);     // get pixel intensity right of current spot
        if (right>0) {               // if part of filament
          endFil = 0;   // filament continues
          newAng = 0; thisTurn = abs(priorAng) - abs(newAng);  // compare to old direction
          if (abs(thisTurn)<abs(turn)) {               // if less of a kink...
            // preLastAng=lastAng; lastAng=newAng; 
            bright=right; nextAng=newAng; turn=thisTurn; nextX=x+1;nextY=y;// save info
          } else if (abs(thisTurn)==abs(turn)) {       // if same kink as before ...
            if (abs(newAng)<abs(nextAng)) {            // only save if closer to horizontal
              // preLastAng=lastAng; lastAng=newAng; 
              bright=right; nextAng=newAng; turn=thisTurn; nextX=x+1;nextY=y; // save info
            } // END IF - more horizontal
          } // END IF - same kink as before
        } // END IF - part of filament

        downright = getPixel(x+1,y+1);   // get pixel intensity below&right of current spot
        if (downright>0) {               // if part of filament
          endFil = 0;   // filament continues
          newAng = 45; thisTurn = abs(priorAng) - abs(newAng);  // compare to old direction
          if (abs(thisTurn)<abs(turn)) {               // if less of a kink...
            // preLastAng=lastAng; lastAng=newAng; 
            bright=downright; nextAng=newAng; turn=thisTurn; nextX=x+1;nextY=y+1;// save info
          } else if (abs(thisTurn)==abs(turn)) {       // if same kink as before ...
            if (abs(newAng)<abs(nextAng)) {            // only save if closer to horizontal
              // preLastAng=lastAng; lastAng=newAng; 
              bright=downright; nextAng=newAng; turn=thisTurn; nextX=x+1;nextY=y+1; // save info
            } // END IF - more horizontal
          } // END IF - same kink as before
        } // END IF - part of filament

        down = getPixel(x,y+1);     // get pixel intensity below current spot
        if (down>0) {               // if part of filament
          endFil = 0;   // filament continues
          newAng = 90; thisTurn = abs(priorAng) - abs(newAng);  // compare to old direction
          if (abs(thisTurn)<abs(turn)) {               // if less of a kink...
            // preLastAng=lastAng; lastAng=newAng; 
            bright=down; nextAng=newAng; turn=thisTurn; nextX=x;nextY=y+1;// save info
          } else if (abs(thisTurn)==abs(turn)) {       // if same kink as before ...
            if (abs(newAng)<abs(nextAng)) {            // only save if closer to horizontal
              // preLastAng=lastAng; lastAng=newAng; 
              bright=down; nextAng=newAng; turn=thisTurn; nextX=x;nextY=y+1; // save info
            } // END IF - more horizontal
          } // END IF - same kink as before
        } // END IF - part of filament

        downleft = getPixel(x-1,y+1);   // get pixel intensity below&left of current spot
        if (downleft>0) {               // if part of filament
          endFil = 0;   // filament continues
          newAng = 135; thisTurn = abs(priorAng) - abs(newAng);  // compare to old direction
          if (abs(thisTurn)<abs(turn)) {               // if less of a kink...
            // preLastAng=lastAng; lastAng=newAng; 
            bright=downleft; nextAng=newAng; turn=thisTurn; nextX=x-1;nextY=y+1;// save info
          } else if (abs(thisTurn)==abs(turn)) {       // if same kink as before ...
            if (abs(newAng)<abs(nextAng)) {            // only save if closer to horizontal
              // preLastAng=lastAng; lastAng=newAng; 
              bright=downleft; nextAng=newAng; turn=thisTurn; nextX=x-1;nextY=y+1; // save info
            } // END IF - more horizontal
          } // END IF - same kink as before
        } // END IF - part of filament

        left = getPixel(x-1,y);     // get pixel intensity left of current spot
        if (left>0) {               // if part of filament
          endFil = 0;   // filament continues
          newAng = 180; thisTurn = abs(priorAng) - abs(newAng);  // compare to old direction
          if (abs(thisTurn)<abs(turn)) {               // if less of a kink...
            // preLastAng=lastAng; lastAng=newAng; 
            bright=left; nextAng=newAng; turn=thisTurn; nextX=x-1;nextY=y;// save info
          } else if (abs(thisTurn)==abs(turn)) {       // if same kink as before ...
            if (abs(newAng)<abs(nextAng)) {            // only save if closer to horizontal
              // preLastAng=lastAng; lastAng=newAng; 
              bright=left; nextAng=newAng; turn=thisTurn; nextX=x-1;nextY=y; // save info
            } // END IF - more horizontal
          } // END IF - same kink as before
        } // END IF - part of filament

        upleft = getPixel(x-1,y-1);   // get pixel intensity above&left of current spot
        if (upleft>0) {               // if part of filament
          endFil = 0;   // filament continues
          newAng = -135; thisTurn = abs(priorAng) - abs(newAng); // compare to old direction
          if (abs(thisTurn)<abs(turn)) {               // if less of a kink...
            // preLastAng=lastAng; lastAng=newAng; 
            bright=upleft; nextAng=newAng; turn=thisTurn; nextX=x-1;nextY=y-1;// save info
          } else if (abs(thisTurn)==abs(turn)) {       // if same kink as before ...
            if (abs(newAng)<abs(nextAng)) {            // only save if closer to horizontal
              // preLastAng=lastAng; lastAng=newAng; 
              bright=upleft; nextAng=newAng; turn=thisTurn; nextX=x-1;nextY=y-1; // save info
            } // END IF - more horizontal
          } // END IF - same kink as before
        } // END IF - part of filament

        if (!endFil) { // filament continues ...
          x = nextX;  y = nextY;               // update position
          preLastAng=lastAng; lastAng=nextAng; // update angles
          // store info
          filX[filLen] = x; filY[filLen] = y;  // store pixel position
          filBright[filLen] = getPixel(x,y);   // store brightness for weight
          filLen++;                            // filament longer
          setPixel(x,y,0);                     // remove pixel to prevent double counting
        } // END IF filament continues 
      } // end while endFil: until end of filament  // // // // // // // // // // // // // // //

      // 3. calculate filemant angles  // // // // // // // // // // // // // // //
      for (i=0; i<filLen; i++) { // all positions along filament
         // calculate angle of filament (smoothing over 'range')
        left = i-range;                       // get array index
        right = i+range;                      // get array index
            // range defines the window over which the angle is calculated.
        if (left<0) { left = 0; }             // fix outlier
        if (right>=filLen) { right=filLen-1; } // fix outlier
        xDiff = filX[right] - filX[left];     // distance in x
        yDiff = filY[right] - filY[left];     // distance in y
        localAngRad = atan2(yDiff,xDiff);     // calculate angle in radians
        if (localAngRad > PI/2) { localAngRad -= PI; }  // correct angles > 90 degrees
        if (localAngRad < -PI/2) { localAngRad += PI; } // correct angles < -90 degrees
        // 4. add up angles for calculation of mean values
        nAng++;                               // one more measurement
        tBright = tBright + filBright[i];     // increase total brightness
        // 5. plot angles in new image
        selectImage(angleRadians);
        setPixel(filX[i], filY[i], localAngRad); // plot angles (in radians) in Angles image 
        selectImage(angleImage);
        setPixel(filX[i], filY[i], localAngRad*127/(PI/2)+128); // plot angles (for 8-bit) in Angles image 
        selectImage(temporary);
      } // end for i: calculate angles for all positions along filament  // // // // // // // // // // // // // // //
    } // END IF ... pixel value > 0 (start of filament)

  } // ENDFOR all y's in image
} // ENDFOR all x's in image

// get raw angle statistics
selectImage(angleRadians); // bring angle image to front
getStatistics(area,angleMean,angleMin,angleMax,angleSD); // get statistics for unweighted angles IN RADIANS

// get weighted angle statistics
imageCalculator("multiply create 32-bit", angleRadians, intensityImage); // new image with angles weighted by brightness
weightedAnglesImage = getImageID();                                    // remember the handle
run("Multiply...", "value="+nAng/tBright);                   // normalize for total brightness, not number of angles
getStatistics(area,angleWmean,angleWmin,angleWmax,angleWsd); // get statistics for weighted angles IN RADIANS

// loop over all pixels to calculate order parameter = deviation from horizontal AND mean angle
orderSum = 0; orderWsum = 0;      // intitialization
for (x=0; x<w; x++) {
  for (y=0; y<h; y++) {
    selectImage(angleRadians);                   // bring 32-bit angle image to front
    localAng = getPixel(x,y);                    // extract local angle 
    localAngRot = localAng - angleWmean;         // ... and subtract weighted mean
    if (!isNaN(localAngRot)) { // if pixel value is not NaN (i.e. if part of filament) ... continue analysis
      order = (2*cos(localAng)*cos(localAng))-1; // calculate order parameter for this angle
      orderRot  = (2*cos(localAngRot)*cos(localAngRot))-1;  // calculate order parameter for rotated angle
      selectImage(intensityImage);      // bring filaments image to front to extract brightness
      brightness = getPixel(x,y);       // get brightness value
      orderW = order * brightness;      // calculate weighted order parameter
      orderWrot = orderRot * brightness;// calculate weighted order parameter for rotated angle
      orderWsum += orderW;               // add up all order parameters
      orderWrotSum += orderWrot;         // add up all rotated order parameters 
    } // END IF ... pixel value not NaN (part of filament)
  } // ENDFOR all y's in image
} // ENDFOR all x's in image
orderWmean = orderWsum/tBright;          // calculate mean weighted order parameter
orderWmeanRot = orderWrotSum/tBright;    // calculate mean weighted order parameter

// calculate mean deviation from horizontal 
selectImage(weightedAnglesImage); 
run("Abs");
getStatistics(area,meanDeviation);   // get statistics for mean deviation from horizontal

selectImage(temporary);
getStatistics(area,mean,min,max);
if (max>0) { // apparently not all actin filaments measured
  showMessage("Alert!","Some actin filaments were not measured!");
} else {
  run("Close");
} // endif: everything okay?


//// 6. save values in parameter arrays

p_count = p_values[0];

p_names[p_count+1] = "Mean Angle (Madison2015)";
p_names[p_count+2] = "Angle Standard Deviation (Madison2015)";
p_names[p_count+3] = "Weighted Mean Angle";
p_names[p_count+4] = "Angular Variation";
p_names[p_count+5] = "Order Parameter (relative to horizontal)";
p_names[p_count+6] = "Order Parameter";
p_names[p_count+7] = "Weighted Angle Deviation";

p_values[p_count+1] = angleMean*180/PI;
p_values[p_count+2] = angleSD*180/PI;
p_values[p_count+3] = angleWmean*180/PI;
p_values[p_count+4] = angleWsd*180/PI;
p_values[p_count+5] = orderWmean;
p_values[p_count+6] = orderWmeanRot;
p_values[p_count+7] = meanDeviation*180/PI;

p_core[0]++; // one more core parameter measured
p_core[p_core[0]] = p_count+3; // Weighted Mean Angle is a core measurement
p_core[0]++; // one more core parameter measured
p_core[p_core[0]] = p_count+4; // Weighted Angle SD is a core measurement
p_core[0]++; // one more core parameter measured
p_core[p_core[0]] = p_count+6; // Weighet Angle Order Parameter (relative to mean weighted angle) is a core measurement

p_values[0] = p_count+7; // Seven new parameters measured


//// 7. cleanup

run("Clear Results");
close("Results");

selectImage(angleRadians); close();
selectImage(weightedAnglesImage); close();

selectImage(angleImage);                
run("Phase");
getLut(reds,greens,blues);
reds[0]=0; greens[0]=0; blues[0]=0;     // set non-filament pixel values to black
setLut(reds,greens,blues);


return(angleImage);

///////////////////////////////////////////////////////////////////////////////

} // END of FUNCTION angles

///////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////

function distances(maskImage, fullFrame, p_names, p_values, p_core) {

///////////////////////////////////////////////////////////////////////////////

// PURPOSE:
// calculates distribution of actin filaments in cells

// INPUT
// maskImage = pointer to image window with skeletonized filaments
// fullFrame = distance calculation in full image (1) or in cell area only (0)
// p_names = 
// p_values
// p_core
// medDist = placeholder for median distance value
// resultsTable = pointer to table that will hold the result values

// RETURN VALUE
// pointer to distance map image
// maximal distance, median distance, skewness of distances
//   (passed back in p_values)
//   (core parameters: median distance)

// PROCEDURE
// 1. invert image, select filaments, create cell outline
// 2. create distance map
// 3. write statistics to arrays


// VERSIONS
// 1.0 = AN.210105
// 1.1 = AN.210108
//     - doing distance map on original skeletonized image
//     - added distribution of filaments (areas, intensities) in cell areas
// 1.2 = AN.210115
//     - added additional table for ratio data
// 1.3 = AN.210218
//     - replaced Left/Right ratios with lower/higher ratios
// 1.4 = AN.210405
//     - added possibility to use predefined cell area
// 2.0 = AN.230207
//     - fixed error with white background in part 2
//     - added analysis of larger "empty" areas: use MaxEntropy threshold for definition
// 2.1 = AN.230808
//     - offer option to calculate distances for entire image, not just cell area (part 2)
// 3.0 = AN.230921
//     - removed cell region analysis
// 3.1 = AN.240610
//     - introduced median distance as pass-back variable
// 4.0 = AN.240613
//     - removed results table and replaced with results arrays

///////////////////////////////////////////////////////////////////////////////

// 1. invert image, select filaments, create cell outline

//// get info of skeletonized image
selectImage(maskImage);
getPixelSize(unit,pxW,pxH);
px = (pxW+pxH)/2; // this is needed since the confocal produces not quite square pixels?
w = getWidth();
h = getHeight();
original = getTitle();
extension = lastIndexOf(original,"-Mask.tif"); // find where extension starts
basename = substring(original,0,extension); // remove extension 

//// create new image for processing
run("Select All"); run("Copy");   // copy image data
newImage("temp","8-bit",w,h,1);   // make new temporary image
setVoxelSize(px,px,1,"micron");   // define pixel size
run("Paste");                     // paste image data into new window
inverted = getImageID();

//// identify filaments, create cell outline
// roiManager("reset");              // forget previous selections
setThreshold(255, 255);           // identify all filaments
run("Create Selection");          // select all filaments
roiManager("add");                // remember this selection
filamentSelection = RoiManager.size - 1; // for easy reference
if (RoiManager.size > 1) { // cell area was predefined
  cellSelection = cellROI;
} else { // cell area not predefined
  run("Convex Hull");               // select minimal area around filaments (= cell)
  roiManager("add");                // remember this selection
  cellSelection = RoiManager.size - 1; // for easy reference
} // ENDIF cell area predefined?
roiManager("deselect");
resetThreshold();


// 2. create distance map

run("Convert to Mask"); // turn image to binary 
run("Invert");
run("Options...", "iterations=1 count=1 black edm=32-bit do=Nothing");
    // this ensures that distances are measured properly
run("Distance Map"); // create new window with distance map
run("Multiply...", "value="+px); // convert pixel distances to micron
setVoxelSize(px,px,1,"micron");   // define pixel size

// calculate maximal and median distance & skewness for entire image frame
getMinAndMax(minDist,maxDistFF);    // maximal distance of pixels from filaments in full image
run("Set Measurements...", "min median skewness redirect=None decimal=3");
run("Measure");
medianDistFF = getValue("Median");  // median distance of pixels from filaments in full image
skewDistFF = getValue("Skew");      // skewness of distance distribution in full image frame

roiManager("select",cellSelection); // reapply cell selection
run("Measure");
maxDistCA = getValue("Max");        // maximal distance of pixels form filaments in cell area
medianDistCA = getValue("Median");  // median distance of pixels from filaments in cell area
skewDistCA = getValue("Skew");      // skewness of distance distribution in cell area

if (!fullFrame) { // IF measure distances in cell area only ===
  setBackgroundColor(0,0,0);
  run("Clear Outside"); // set everything outside of cell to zero
} // END IF measure distances in cell area only ===
setMinAndMax(0,10); // enhance contrast in distance map image (from 0 to 10 microns)
run("Fire"); // apply LUT "Fire"
rename(basename+"-DistanceMap.tif");
distMap = getImageID(); // remember distance map image


// 3. write results into arrays.

p_count = p_values[0]; // number of paramters already measured

p_names[p_count+1] = "Maximal Distance";
p_names[p_count+2] = "Median Distance";
p_names[p_count+3] = "Skewness of Distances";

if (fullFrame) { // IF Distance measurements based on full frame **
  p_values[p_count+1] = maxDistFF;
  p_values[p_count+2] = medianDistFF;
  p_values[p_count+3] = skewDistFF;
} else { // ELSE Distance measurements based on cell area **
  p_values[p_count+1] = maxDistCA;
  p_values[p_count+2] = medianDistCA;
  p_values[p_count+3] = skewDistCA;
} // END IF full Frame or cell area?  **

p_values[0] = p_count + 3; // three more measurements

p_core[0]++;
p_core[p_core[0]] = p_count+2; // median distance is core parameter


// cleanup

run("Clear Results");
close("Results");

roiManager("reset");              // forget selections 
close("ROI Manager");
close("temp");                    // close temporary window

return distMap; // return handle of density map image


///////////////////////////////////////////////////////////////////////////////

} // END of FUNCTION distances

///////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////

function table_stats(theTable, theTablePlusStats) {

///////////////////////////////////////////////////////////////////////////////

// PURPOSE
// Calculate standard statistics for all columns in a table:
//   mean, sd, CoV, median, min, max, range, relative range

// INPUT:
// theTable = pointer to table with numerical values
// theTablePlusStats = pointer to new table with statistics

// RETURN VALUE
// number of processed files

// PROCEDURE
// 1. get all column headings and find text columns
// 2. turn column into array
// 3. do calculations and write at bottom of column
// 4. save everything (data + stats) to new column

// VERSIONS:
// 1.0 = AN.210115

///////////////////////////////////////////////////////////////////////////////

// 1. get all column headings and find text columns

headers = split(Table.headings(theTable),"\t");
// num = newArray(headers.length);
stringColumns = 0;
for (i=0; i< headers.length; i++) { // all columns **

// 2. turn column into array
  data = Table.getColumn(headers[i],theTable); // original data
  newData = newArray(data.length+9); // new array for data + statistics
  for (j=0; j<data.length; j++) { newData[j] = data[j]; } // copy data
  if (isNaN(parseFloat(data[0]))) { // check if column contains numbers ===
    stringColumns++;
    newData[data.length] = "xxEMPTYxx"; // empty row between data and stats
    newData[data.length+1] = "Mean";
    newData[data.length+2] = "StdDev";
    newData[data.length+3] = "Coefficient of Variance";
    newData[data.length+4] = "Median";
    newData[data.length+5] = "Min";
    newData[data.length+6] = "Max";
    newData[data.length+7] = "Range";
    newData[data.length+8] = "Relative Range";
  } else {  // numerical column ===

// 3. do calculations and write at bottom of column
    Array.getStatistics(data, min, max, mean, stdDev);
    Array.sort(data);
    middle = (data.length+1)/2;
    median = (data[floor(middle)-1] + data[Math.ceil(middle)-1]) / 2; 
    cov = stdDev/mean;
    range = max-min;
    relRange = range/median;

    newData[data.length+1] = mean;
    newData[data.length+2] = stdDev;
    newData[data.length+3] = cov;
    newData[data.length+4] = median;
    newData[data.length+5] = min;
    newData[data.length+6] = max;
    newData[data.length+7] = range;
    newData[data.length+8] = relRange;
  } // END ELSE numerical column ===
  Table.setColumn(headers[i],newData,theTablePlusStats);
  Table.set(headers[i],data.length,"",theTablePlusStats); // clear empty row

} // END FOR all columns **

// 4. finish

return data.length;

///////////////////////////////////////////////////////////////////////////////

} // END FUNCTION table_stats

///////////////////////////////////////////////////////////////////////////////
