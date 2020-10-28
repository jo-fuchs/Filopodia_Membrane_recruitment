// Plasma membrane recruitment with combined Filopodia analysis
//
// Filopodia-analysis: Counts Filopodia around a cell
// 		as intensity maxima on a line 1 µm outside of cell mask
// 
//	PM-recruitment-analysis
//  	Measures intensity at cell surface, intracellular & perinuclear
//		creates 1µm thick band to measure intensity (mean & median) in all regions
//    
// requires:
//		membrane marker (fGFP or F-actin) optionally with Dapi (for perinuclear intensity)
//		manual selection of cell of interest in ROI-manager
//
// Output:
//    ROIs saved as Image_ROIname_ROIs.zip
//    Intensities saved as Image_ROIname_PMrecruitment.csv
//	  Filopodia density saved as Image_ROIname_Filopodia.csv
//
// v0.1 26.10.2020 Joachim Fuchs
//



// define channels
GFPchannel = 1;
PRG2channel = 4;
Dapichannel = 3; // if no Dapi imaged, set to 0




// some general setup
run("Clear Results");
run("Input/Output...", "jpeg=100 gif=-1 file=.csv save_column");
run("Set Measurements...", "area mean standard median display redirect=None decimal=2");

path = getDirectory("image");
dirimg = path + "Measurements/";
File.makeDirectory(dirimg); 
imgName = getTitle();
run("Line Width...", "line="+5); 
run("Select None");
roiManager("deselect");


// refine cell ROI
run("Duplicate...", "duplicate");
roiManager("Select", 0);
rName = Roi.getName;


Stack.setChannel(GFPchannel);
run("Median...", "radius=2 slice");
setAutoThreshold("Huang dark");
run("Analyze Particles...", "size=200-Infinity clear include add slice");
close();
roiManager("Select", 0);
Stack.setChannel(PRG2channel);
cell_surface = getValue("Area");
mean_intensity = getValue("Mean");
Stack.setChannel(GFPchannel);

// remove some long thin processes, enlarge by 1 µm
run("Enlarge...", "enlarge=-1");
run("Enlarge...", "enlarge=2");
	// remove unconnected parts & fill holes
	run("Create Mask");
	run("Analyze Particles...", "size=200-Infinity include add");
	run("Close");
	roiManager("Select", 1);

run("Area to Line");
roiManager("add");
roiManager("Select", 2);
roiManager( "Rename", "Filopodia" )
Length = getValue("Length");

/// have an adaptable threshold depending on intensities of the lines
roiManager("measure");
selectWindow(imgName);
intens = getResult("Median", 0);


run("Plot Profile");
run("Find Peaks", "min._peak_amplitude=" + intens + 
	" min._peak_distance=1.5 min._value=NaN max._value=0 exclude list");

//Get number of maxima
selectWindow("Plot Values");
a = Table.getColumn("X1");
Filos = lengthOf(a);
run("Close");

dens = Filos / Length;

// Create Results
run("Clear Results");
setResult("Image", 0, imgName);
setResult("Channel", 0, GFPchannel);
setResult("Area", 0, cell_surface);
setResult("MeanIntensity", 0, mean_intensity);
setResult("Circumference", 0, Length);
setResult("FiloNumber", 0, Filos); 
setResult("FiloDensity", 0, dens); 

saveAs("Results", dirimg + imgName + "_" + rName + "_Filopodia.csv");


// PM-localization

selectWindow(imgName);
Stack.getPosition(channel, slice, frame)
run("Select None");
roiManager("deselect");






// Plasma membrane Selection
roiManager("Select", 0);
Stack.setPosition(GFPchannel, slice, frame);
run("Enlarge...", "enlarge=-1.2");
run("Make Band...", "band=1");
roiManager("add");
roiManager("Deselect");


// Intracellular selection
run("Select None");
roiManager("Select", 0);
run("Enlarge...", "enlarge=-3");
	// Secure against unconnected ROIs remove unconnected parts & fill holes
	run("Create Mask");
	run("Analyze Particles...", "size=100-Infinity include add");
	run("Close");
	
roiManager("Select", 4);
run("Make Band...", "band=1");
roiManager("add");

roiManager("Select", 4);
roiManager("delete");


if(Dapichannel != 0) {
	// perinuclear selection
	run("Select None");
	run("Duplicate...", "duplicate");
	roiManager("Select", 0);
	Stack.setPosition(Dapichannel, slice, frame)  // Dapi channel
	run("Median...", "radius=2 slice");
	setAutoThreshold("Moments dark");
	run("Analyze Particles...", "size=70-Infinity add include slice");
	
	// This fails for multiple nuclei > merge multiple nuclei before
	n = roiManager("Count");
	a = Array.getSequence(n);
	a = Array.slice(a,5,n);
	if(a.length > 1) { // only if multiple nuclei
		roiManager("select", a);
		roiManager("Combine");
		run("Convex Hull");
		roiManager("add");
		
		//remove individual ROIs
		roiManager("select", a);
		roiManager("delete");
	}
	
	// create band around nuclei
	roiManager("Select", 5);
	run("Enlarge...", "enlarge=-1");
	run("Convex Hull");
	run("Make Band...", "band=1");
	roiManager("add");
	// remove nucleus area
	roiManager("Select", 5);
	roiManager("delete");
	close(); // duplicated image
}

// Measure all
roiManager("Select", 3);
roiManager( "Rename", "Membrane" )
Stack.setPosition(GFPchannel, slice, frame)
run("Measure");
Stack.setPosition(PRG2channel, slice, frame)
run("Measure");

roiManager("Select", 4);
roiManager( "Rename", "Intracellular" )
Stack.setPosition(GFPchannel, slice, frame)
run("Measure");
Stack.setPosition(PRG2channel, slice, frame)
run("Measure");

if(Dapichannel != 0) {
	roiManager("Select", 5);
	roiManager( "Rename", "Perinuclear" )
	Stack.setPosition(GFPchannel, slice, frame)
	run("Measure");
	Stack.setPosition(PRG2channel, slice, frame)
	run("Measure");
}

// save results
roiManager("select", 1);
roiManager("delete");
roiManager("Deselect");
roiManager("Save", dirimg + imgName + "_" + rName +"_ROIs.zip");
saveAs("Results", dirimg + imgName + "_" + rName + "_PMrecruitment.csv");


// clean up
run("Close All");

list = getList("window.titles");
 for (i=0; i<list.length; i++){
 winame = list[i];
  selectWindow(winame);
 run("Close");
 }