run("Split Channels");

//get image name
selectImage(1);
var imageName1 = getTitle();
var imageName2 = substring(imageName1, 3, lengthOf(imageName1)-8);


// select 647
selectImage(1);
run("Duplicate...", "title=647");
run("Subtract Background...", "rolling=50");
// RS FISH: max 20,000, sigma 1.50, threshold 0.0020
run("RS-FISH", "image=647 mode=Advanced anisotropy=1.0000 robust_fitting=[No RANSAC] add image_min=0 image_max=20000 sigma=1.5 threshold=0.0014 support=3 min_inlier_ratio=0.10 max_error=1.50 spot_intensity_threshold=289.01 background=[No background subtraction] background_subtraction_max_error=0.05 background_subtraction_min_inlier_ratio=0.10 results_file=[]");
var count = roiManager("count");
roiManager("select", Array.getSequence(count));
roiManager("Combine");
roiManager("Add");
roiManager("Delete");
roiManager("Deselect");
roiManager("Select", 0);
roiManager("Rename", "647");
roiManager("Save", "Y:/People/current/Giacomo/HCR/chir experiment/time course experiment 3/axin and fgf10/reprocessed/rois/647-"+imageName2+".roi");
roiManager("delete");
selectWindow("smFISH localizations");
save("Y:/People/current/Giacomo/HCR/chir experiment/time course experiment 3/axin and fgf10/reprocessed/smFISH Localizations/647-"+imageName2+"-smFISH localizations.txt");
run("Close");
selectWindow("647");
saveAs("tif", "Y:/People/current/Giacomo/HCR/chir experiment/time course experiment 3/axin and fgf10/reprocessed/tifs/647-"+imageName2+".tif");



// select 546
selectImage(2);
run("Duplicate...", "title=546");
run("Subtract Background...", "rolling=50");
// RS FISH: max 20,000, sigma 1.50, threshold 0.0020
run("RS-FISH", "image=546 mode=Advanced anisotropy=1.0000 robust_fitting=[No RANSAC] add image_min=0 image_max=20000 sigma=1.5 threshold=0.0018 support=3 min_inlier_ratio=0.10 max_error=1.50 spot_intensity_threshold=289.01 background=[No background subtraction] background_subtraction_max_error=0.05 background_subtraction_min_inlier_ratio=0.10 results_file=[]");
var count = roiManager("count");
roiManager("select", Array.getSequence(count));
roiManager("Combine");
roiManager("Add");
roiManager("Delete");
roiManager("Deselect");
roiManager("Select", 0);
roiManager("Rename", "546");
roiManager("Save", "Y:/People/current/Giacomo/HCR/chir experiment/time course experiment 3/axin and fgf10/reprocessed/rois/546-"+imageName2+".roi");
roiManager("delete");
selectWindow("smFISH localizations");
save("Y:/People/current/Giacomo/HCR/chir experiment/time course experiment 3/axin and fgf10/reprocessed/smFISH Localizations/546-"+imageName2+"-smFISH localizations.txt");
run("Close");
selectWindow("546");
saveAs("tif", "Y:/People/current/Giacomo/HCR/chir experiment/time course experiment 3/axin and fgf10/reprocessed/tifs/546-"+imageName2+".tif");



// deep cell

selectImage(4);
run("Duplicate...", "title=dapi");
// get dimensions of image
var width = getWidth();
var height = getHeight();
run("Bin...", "x=10 y=10 bin=Min");

selectImage(3);
run("Duplicate...", "title=green");
run("Normalize Local Contrast", "block_radius_x=40 block_radius_y=40 standard_deviations=2 center stretch");
run("Bin...", "x=10 y=10 bin=Median");

run("Merge Channels...", "c1=dapi c2=green create keep");

run("Submit Active Image", "deepcell=http://deepcell.org status=5000 job=3600 select=mesmer");
setOption("ScaleConversions", true);
run("16-bit");
run("glasbey inverted");
run("Scale...", "x=- y=- width=["+width+"] height=["+height+"] interpolation=None average create title=DeepCellOutput.tif");
save("Y:/People/current/Giacomo/HCR/chir experiment/time course experiment 3/lef1 and fgf8/reprocessed/deep cell/"+imageName1+"-deep cell.tif");

selectWindow("Composite.tif");
run("Scale...", "x=- y=- width=["+width+"] height=["+height+"] interpolation=None average create title=DapiGreen.tif");
save("Y:/People/current/Giacomo/HCR/chir experiment/time course experiment 3/lef1 and fgf8/reprocessed/DapiGreen/"+imageName1+"-DapiGreen.tif");

run("Close All");
