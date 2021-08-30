/*
 * - - - MakeIsotropicM40images.ijm - - - 
 * 
 * Author:	Keisuke Ishihara
 * Date:	May 13, 2021
 * 
 * Macro for batch-processing .vsi image files from Olympus spinning disc microscopes for downstream 3D analysis.
 * Output is isotropic z-stacks. Each channel results in a different tif file.
 * CLIJ's makeIsotropic function is used to quickly downsample images on the GPU.
 * 
 * Usage:
 * 1. Prepare output folder.
 * 2. Run the macro on ImageJ/Fiji and follow the instructions in the window.
 * 
 */

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".vsi") suffix
#@ Integer (label = "Nr. of files to process (-1 to process all).", min=-1, value = 2) n2process

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);

	print("------------------------");

	run("CLIJ2 Macro Extensions", "cl_device=");
	Ext.CLIJ2_getGPUProperties(GPU_name, global_memory_in_bytes, OpenCL_version);
	print("using "+ GPU_name + " GPU with " + (global_memory_in_bytes / 1024 / 1024 / 1024)+"GB Memory, OpenCL version " + OpenCL_version);
	
	print("Input directory: " + input );
	print("Output directory: " + output);

	// count number of files ending with suffix
	n_suffix = 0;
	for (i = 0; i < list.length; i++) {
		if(endsWith(list[i], suffix)) {
			n_suffix += 1;
		}
	}
	print(toString(n_suffix)+" files ending with "+ suffix);

	if (n2process>n_suffix) {
		n2process = n_suffix;
	}
	if (n2process<0) {
		n2process = n_suffix;
	}

	print("Process "+toString(n2process)+" files.");

	start = getTime();
	cc = 0;
	n_processed = 0;
	for (i = 0; i < list.length; i++) {
		if(endsWith(list[i], suffix)) {
			cc += 1;
			if (cc <= n2process) {
				print("File " + toString(cc) + " of " + toString(n_suffix) +": "+ list[i]);
				processFile(input, output, list[i]);
				n_processed += 1;
			}
		}
	}
	print("Processed "+toString(n_processed)+" files.");
	print("Execution time: " + toString((getTime()-start)/1000) + " seconds.");  
}

function processFile(input, output, file) {
	// print("Input: " + input + File.separator + file);
	// print("Output: " + output);
	
	// open image as Bio-Formats
	run("Bio-Formats", "open=["+input+File.separator+file+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
	rename("stack");
	getDimensions( width, height, channels, slices, frames );

	run("Split Channels");
	
	// make isotropic
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
		Ext.CLIJ_convertUInt16("isotemp", "isotemp1"); // for some reason this needs to CLIJ not CLIJ2, CLIJx
		Ext.CLIJ2_pull("isotemp1");
		
		saveAs("Tiff", output + File.separator + replace(file, suffix, "") + "_isotropic_C"+ toString(j) +".tif");
		
		Ext.CLIJ2_clear();
	}
	run("Close All");
}
