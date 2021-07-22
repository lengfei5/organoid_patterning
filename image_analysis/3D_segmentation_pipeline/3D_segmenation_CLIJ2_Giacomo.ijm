setBatchMode(false);

// Split channels
run("Split Channels");

selectImage(1);
run("Subtract...", "value=680 stack");
rename("647");

selectImage(2);
image_1 = getTitle();
run("Duplicate...", "title=546 duplicate");
run("Subtract...", "value=624 stack");

// Init GPU
run("CLIJ2 Macro Extensions", "cl_device=");
Ext.CLIJ2_clear();

// make mask for 546 channel
Ext.CLIJ2_push(image_1);

// Gaussian Blur3D
sigma_x = 3.0;
sigma_y = 3.0;
sigma_z = 3.0;
Ext.CLIJ2_gaussianBlur3D(image_1, image_2, sigma_x, sigma_y, sigma_z);
Ext.CLIJ2_release(image_1);
Ext.CLIJ2_pull(image_2);
Ext.CLIJ2_release(image_2);
rename("fgf8_blurred");
Ext.CLIJ2_clear();

// make mask for 488 channel
selectImage(3);
image_1 = getTitle();
Ext.CLIJ2_push(image_1);

// Gaussian Blur3D
sigma_x = 3.0;
sigma_y = 3.0;
sigma_z = 3.0;
Ext.CLIJ2_gaussianBlur3D(image_1, image_2, sigma_x, sigma_y, sigma_z);
Ext.CLIJ2_release(image_1);
Ext.CLIJ2_pull(image_2);
Ext.CLIJ2_release(image_2);
rename("gfp_blurred");
Ext.CLIJ2_clear();


imagePath = "/groups/tanaka/People/current/Giacomo/HCR/chir experiment/time course experiment 4/masks/";
imageName = substring(image_1, 3, lengthOf(image_1)-4);
// thresholding to make masks
selectWindow("gfp_blurred");
setThreshold(400.0800, 1000000000000000000000000000000.0000);
run("Convert to Mask", "method=Default background=Dark black");
//saveAs("tiff", imagePath+"blood-"+imageName);
rename("blood mask");


selectWindow("fgf8_blurred");
setThreshold(1000.0800, 1000000000000000000000000000000.0000);
run("Convert to Mask", "method=Default background=Dark black");
//saveAs("tiff", imagePath+"fgf8-"+imageName);
rename("fgf8 mask");


// subtract the blood
imageCalculator("Subtract create stack","fgf8 mask","blood mask");
selectWindow("Result of fgf8 mask");
rename("Fgf8 minus blood");


// add masks to 3D manager
selectWindow("Fgf8 minus blood");
run("3D Manager");
Ext.Manager3D_AddImage();

imagePath = "/groups/tanaka/People/current/Giacomo/HCR/chir experiment/time course experiment 4/results/";

// make measurements and save them for 546
Ext.Manager3D_Select(0);
Ext.Manager3D_DeselectAll();
selectWindow("546");
Ext.Manager3D_Select(0);
Ext.Manager3D_Quantif();
Ext.Manager3D_SaveResult("Q",imagePath+"546-"+imageName+".csv");
Ext.Manager3D_CloseResult("Q");
Ext.Manager3D_List();
Ext.Manager3D_SaveResult("L",imagePath+"546-"+imageName+".csv");
Ext.Manager3D_CloseResult("L");

// make measurements and save them for 647
Ext.Manager3D_Select(0);
Ext.Manager3D_DeselectAll();
selectWindow("647");
Ext.Manager3D_Select(0);
Ext.Manager3D_Quantif();
Ext.Manager3D_SaveResult("Q",imagePath+"647-"+imageName+".csv");
Ext.Manager3D_CloseResult("Q");
Ext.Manager3D_List();
Ext.Manager3D_SaveResult("L",imagePath+"647-"+imageName+".csv");
Ext.Manager3D_CloseResult("L");

// close 3d roi manager
Ext.Manager3D_Delete();
Ext.Manager3D_Close();

run("Close All");
