/*
 * Macro to process multiple images (Olympus .vsi) in a folder
 * uses CLIJ to bin the images for isotropic
 */

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".vsi") suffix
//#@ Integer (label = "Nr. of files to process.", value = 1) n2process

// See also Process_Folder.py for a version of this code
// in the Python scripting language.

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);

	print("Input directory: " + input );
	print("Output directory: " + output);

	// count number of files ending with suffix
	n_suffix = 0;
	for (i = 0; i < list.length; i++) {
		if(endsWith(list[i], suffix)) {
			n_suffix += 1;
		}
	}

	n_processed = 0;
//	for (i = 0; i < list.length; i++) {
	for (i = 0; i < 5; i++) {
		if(endsWith(list[i], suffix)) {
			n_processed += 1;
//			print("File " + toString(i) + " of " + toString(list.length) +": "+ list[i]);
			print("File " + toString(n_processed) + " of " + toString(n_suffix) +": "+ list[i]);
			processFile(input, output, list[i]);
		}
	}
	print("Processed"+toString(n_processed)+" files.");
}

function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	// print("Input: " + input + File.separator + file);
	// print("Output: " + output);
	
	// open image as Bio-Formats
	run("Bio-Formats", "open=["+input+File.separator+file+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
	rename("stack");
	getDimensions( width, height, channels, slices, frames );

	run("Split Channels");
	
	// make isotropic
	run("CLIJ2 Macro Extensions", "cl_device=[Quadro P1000]");

	for (j = 1; j < (channels+1); j++) {
		input_image = "C"+toString(j)+"-stack";
		Ext.CLIJ2_push(input_image);
		isotropic_image = "isotemp";
		original_voxel_size_x = 0.65;
		original_voxel_size_y = 0.65;
		original_voxel_size_z = 3.0;
		new_voxel_size = 3.0;
		Ext.CLIJx_makeIsotropic(input_image, isotropic_image, original_voxel_size_x, original_voxel_size_y, original_voxel_size_z, new_voxel_size);
		
		// convert 32 to 16 bits
		//intensityScalingFactor = 1.0;
		//Ext.CLIJ2_multiplyImageAndScalar("isotemp", "isotemp1", intensityScalingFactor);
		//Ext.CLIJ2_convertUInt16("isotemp1", "isotemp2");
		//Ext.CLIJ2_pull("isotemp2");
		Ext.CLIJ_convertUInt16("isotemp", "isotemp1");
		Ext.CLIJ2_pull("isotemp1");
		
		saveAs("Tiff", output + File.separator + replace(file, suffix, "") + "_isotropic_C"+ toString(j) +".tif");
		Ext.CLIJ2_clear();
	}
	run("Close All");
}
