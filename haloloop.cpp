#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2\opencv.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <direct.h>
#include <ctime>




using namespace cv;
using namespace std;

//Purpose: analyse images taken by allsky camera for halo phenomena
//authors: Stephen Sorenson (6-2013 to 4-2014), Shelby Richard (5-2014 to 2015), Michelle King (2015-17), Morton Greenslit (2015-27), Sylke Boyd (6-2013 to present)
//version: based on the 2015 version of haloloop for UMM allsky images
// 06-12-2018:  reworking for ARM TSI images


int main(int argc, char** argv)
{
	clock_t startTime = clock();//keeping track of runtime
	stringstream sc;

	
	//header to screen
	cout << " ************************ Program haloloop *********************** \n";
	cout << "     Purpose: search ice halos in series of allsky images \n";
	cout << " authors: S. Sorenson, S. Richard, M. King, M. Greenslit, S.Boyd   \n";
	cout << "                University of Minnesota - Morris \n";
	cout << "    version 2018 - processing TSI images from NOAA ARM sites       \n";
	cout << " ***************************************************************** \n";
	cout << "\n";
	
	//**************  Reading the paths and names to the images and creating processfile ************************
		//open a file with the list of folders and imagefile names
		cout << "* Looking for the text file that contains the list of images to be processed. \n";
		cout << "*     This file can be produced in a linux shell using the command: \n";
		cout << "*     ls -R1  foldernames  > imagelist \n";
		cout << "* File name: ";
		string imagelistfile;
		cin >> imagelistfile;
		
		//check for existence of imagelist 
		const char * c = imagelistfile.c_str();
		if (access( c, R_OK ) != 0) 
		{ 
			cout << "!! File " << imagelistfile << " not found! Program haloloop stops here !!\n"; exit(0); 
		} else { cout << "*   Found file " << imagelistfile << " \n";}
		
		//open the imagelist and prepare the temporary processlist
		ifstream imagelist(imagelistfile);
		ofstream processlistout("haloloop-processlist.txt"); //a temporary file to contain the cleaned-up paths and names of the images for the carpet
		
		//open a log file that will contain all parameters for this run
		string runlogfile;
		sc << imagelistfile << "-haloloop-out-RunLog.txt"; runlogfile = sc.str(); sc.str("");//assign the name of the parameterfile
		ofstream runlog(runlogfile); //messages and running output
		runlog << "********* " << runlogfile << " ******************\n";
		runlog << "**  will hold confirmation of parameters, usable as parameter file\n";
		cout << "*   all progress recorded in " << runlogfile << "  \n";
		
		//read the folder and filenames, identify origin and time, and store them in a temporary processlist
		string CurrentImage, CurrentFolder, dummystring;
		int imagecounter = 0;//total image counter
		while (imagelist){
			getline(imagelist,CurrentFolder); //reads the first current folder name
			CurrentFolder = CurrentFolder.erase(CurrentFolder.size()-1); //erases the colon behind the string
			int anotherfile = 1; //decides whether all files in a folder have been read
			int imagesinfolder = 0; //counts the images in this folder
			while (anotherfile == 1){
				getline(imagelist,dummystring);
				if (dummystring != "") //another filename has been read (not empty)
				{
					CurrentImage = dummystring;
					string originstr = CurrentImage.substr(0,3);
					string yearstr =   CurrentImage.substr(20,4);
					string monthstr =  CurrentImage.substr(24,2);
					string daystr =    CurrentImage.substr(26,2);
					string hourstr =   CurrentImage.substr(29,2);
					string minutestr = CurrentImage.substr(31,2);
					string secondstr = CurrentImage.substr(33,2);
					//writing information out to the processlist
					processlistout << CurrentFolder << "/" << CurrentImage << "\t" << originstr << "\t" << yearstr << "\t" << monthstr << "\t" << daystr << "\t" << hourstr << "\t" << minutestr << "\t" << secondstr << "\n";
					imagecounter++;//advance the counters
					imagesinfolder++;
				}	
				else { anotherfile = 0;} //an empty line was encountered at the end of the list for one folder
			} //end of imagelist for a single folder
			//some information to standard output
			cout << "*   Haloloop will process folder \t" << CurrentFolder << "\t with \t" << imagesinfolder << "\t images \n" ;
		} //the whole imagelist has been read
		// define total number of images to be processed
		int totalimagecount = imagecounter;
		runlog << "**   Haloloop will process a total of " << totalimagecount << " images. \n";
		imagelist.close();
		processlistout.close(); //close this to open it again as an ifstream
	//*************done Reading the paths and names to the images and creating processfile ************************
	
	
	//***************** Reading the input parameter file ************************************************************
		double Pi = 3.14159265359;		
		//open the input-parameter file
			string parameterfilename;
			sc << imagelistfile << "-haloloop-in-prms.txt"; parameterfilename = sc.str(); sc.str("");//assign the name of the parameterfile
		//check for existence of parameterfile
			c = parameterfilename.c_str();
			if (access( c, R_OK ) != 0) 
			{ 
				cout << "!! File " << parameterfilename  << " not found! Program haloloop stops here !!\n"; exit(0); 
			} else { cout << "*   Found file " << parameterfilename << " \n";}
		
			ifstream parameterfile(parameterfilename);	
			for (int iline=0; iline <3; iline++) {getline(parameterfile,dummystring);} //first three lines are comments
		//Image pre-processing: 
			getline(parameterfile,dummystring);//skip a comment line used as header in parameterfile
			runlog << dummystring << " \n";
			double RelcolShift[4]; //relative color shift corrections, [1] B shift, [2] G shift, [3] R shift, [4] alpha
				parameterfile >> RelcolShift[0] >> RelcolShift[1] >> RelcolShift[2] >> RelcolShift[3] ; getline(parameterfile,dummystring);
				runlog << "\t" << RelcolShift[0] << "\t" << RelcolShift[1] << "\t" <<  RelcolShift[2] << "\t" <<  RelcolShift[3] << dummystring << " \n"; 			
			Point2f HCenterShift; //corrections for the horizon ellipse: centerposition; These may change every time an adjustment to the TSI was made.
				parameterfile >> HCenterShift.x >> HCenterShift.y; getline(parameterfile,dummystring); //the getline ends the line in the parameterfile while ignoring it(allows comments in parameterfile)
				runlog << "\t" << HCenterShift.x << "\t"<< HCenterShift.y << dummystring << " \n"; //this format allows that the same line could be copied back into the parameterfile
			double NorthAlignAngle; //angle alignment corrections
				parameterfile >> NorthAlignAngle; getline(parameterfile,dummystring);
				runlog << "\t" << NorthAlignAngle << dummystring << " \n"; 			
			double HXadjust; //corrections for the horizon ellipse: scaling of horizontal axis relative to image width
				parameterfile >> HXadjust; getline(parameterfile,dummystring);
				runlog << "\t" << HXadjust << dummystring << " \n"; 
			double eps; //corrections for the horizon ellipse: ellipticiy  = yaxis/xaxis;
				parameterfile >> eps; getline(parameterfile,dummystring);
				runlog << "\t" << eps << dummystring << " \n"; 
			double MaximumZenithAngle; //defines the maximum zenith angle available for analysis in the images. This may vary by location		
				parameterfile >> MaximumZenithAngle; getline(parameterfile,dummystring);
				runlog << "\t" << MaximumZenithAngle << dummystring << " \n"; 
			int SaveHorizonMat; //save original image with marked original horizon for every image (avoid that for big runs!)
				parameterfile >> SaveHorizonMat; getline(parameterfile,dummystring);
				if (totalimagecount > 200 ){SaveHorizonMat =0;} //don't allow excess images for large runs - avoid the pain
				runlog << "\t" << SaveHorizonMat << dummystring << " \n"; 
			int SaveCircImageMat; //if set to 1 saving the image as stretched to circular horizon as circimage (avoid that for big runs!)
				parameterfile >> SaveCircImageMat; getline(parameterfile,dummystring);
				if (totalimagecount > 200 ){SaveCircImageMat =0;} //don't allow excess images for large runs - avoid the pain
				runlog << "\t" << SaveCircImageMat << dummystring << " \n"; 
			int SaveNorthAlignMat; //if set to 1 saving the circular image as rotated to north (avoid that for big runs!)
				parameterfile >> SaveNorthAlignMat; getline(parameterfile,dummystring);
				if (totalimagecount > 200 ){SaveNorthAlignMat =0;} //don't allow excess images for large runs - avoid the pain
				runlog << "\t" << SaveNorthAlignMat << dummystring << " \n"; 
			int SaveSunMarked; //if set to 1 saving the circle-corrected image with sun and center marked is saved (avoid that for big runs!) 
				parameterfile >> SaveSunMarked; getline(parameterfile,dummystring);
				if (totalimagecount > 200 ){SaveSunMarked =0;} //don't allow excess images for large runs - avoid the pain
				runlog << "\t" << SaveSunMarked << dummystring << " \n"; 
			int SaveArcImageMat; //if set to 1 saving the image with arc projection (avoid that for big runs!)
				parameterfile >> SaveArcImageMat; getline(parameterfile,dummystring);
				if (totalimagecount > 200 ){SaveArcImageMat =0;} //don't allow excess images for large runs - avoid the pain
				runlog << "\t" << SaveArcImageMat << dummystring << " \n"; 
			int SaveSkyMaskMat; //if set to 1 the final prepared image with shadows masked is saved- this one is the basis for the local sky map (avoid that for big runs!) 	
				parameterfile >> SaveSkyMaskMat; getline(parameterfile,dummystring);
				if (totalimagecount > 200 ){SaveSkyMaskMat =0;} //don't allow excess images for large runs - avoid the pain
				runlog << "\t" << SaveSkyMaskMat << dummystring << " \n"; 
			int SaveLocalSkyMat; //if set to 1 the local sky map surrounding the sun is saved (avoid that for big runs!) 	
				parameterfile >> SaveLocalSkyMat; getline(parameterfile,dummystring);
				if (totalimagecount > 500 ){SaveLocalSkyMat =0;} //don't allow excess images for large runs - avoid the pain
				runlog << "\t" << SaveLocalSkyMat << dummystring << " \n"; 
			double ShadowStripWidth, ShadowStripOffset; // the width in pixels used to mask the shadow strip; in some locations/times an offset seems necessary to center the shadowstrip mask.
				parameterfile >> ShadowStripWidth >> ShadowStripOffset; getline(parameterfile,dummystring);
				runlog << "\t" << ShadowStripWidth << "\t" << ShadowStripOffset << dummystring << " \n"; 
			int LSMCropSizesd; // 	50				LSMCropSizesd  degree side length of Local Sky Mat, 50 is ok (halo diameter 44 deg, plus edge)
				parameterfile >> LSMCropSizesd; getline(parameterfile,dummystring);
				runlog << "\t" << LSMCropSizesd << dummystring << " \n"; 
		//Radial Intensity Analysis
			getline(parameterfile,dummystring);//skip a comment line used as header in parameterfile
			runlog << dummystring << " \n";
			int SolDistBin; //bin width centered at each analysis radius to build the average (overlap)
				parameterfile >> SolDistBin; getline(parameterfile,dummystring);
				runlog << "\t" << SolDistBin << dummystring << " \n"; 
			int SaveAveRadInt; //if 1 the data table for the average radial intensities (ARI) for the current image will be saved (avoid that for big runs!) 
				parameterfile >> SaveAveRadInt; getline(parameterfile,dummystring);
				if (totalimagecount > 200 ){SaveAveRadInt =0;} //don't allow excess images for large runs - avoid the pain
				runlog << "\t" << SaveAveRadInt << dummystring << " \n"; 
		//Halo search using running average analysis
			getline(parameterfile,dummystring);//skip a comment line used as header in parameterfile
			runlog << dummystring << " \n";
			int RALowLim; //lowest radial distance for which RA and its derivative will be computed
				parameterfile >> RALowLim; getline(parameterfile,dummystring);
				runlog << "\t" << RALowLim << dummystring << " \n"; 
			int RAHighLim; //highest radial distance for which RA and its derivative will be computed
				parameterfile >> RAHighLim; getline(parameterfile,dummystring);
				runlog << "\t" << RAHighLim << dummystring << " \n"; 
		//line fit and line fit area analysis
			getline(parameterfile,dummystring);//skip a comment line used as header in parameterfile
			runlog << dummystring << " \n";
			double rfitlimmin;  // low radial boundary for linear fit in skydegrees, 15 degrees good
				parameterfile >> rfitlimmin; getline(parameterfile,dummystring);
				runlog << "\t" << rfitlimmin << dummystring << " \n"; 
			double rfitlimmax;  // high radial boundary for linear fit in skydegrees, 29 good
				parameterfile >> rfitlimmax; getline(parameterfile,dummystring);
				runlog << "\t" << rfitlimmax << dummystring << " \n"; 
			int SaveLineFitPRMs; //creates file with a line for each processed image, containing sloe and intercept for each quadrant (can be done for big runs)
				parameterfile >> SaveLineFitPRMs; getline(parameterfile,dummystring);
				runlog << "\t" << SaveLineFitPRMs << dummystring << " \n"; 
		// Handling SkyTypeScores
			getline(parameterfile,dummystring);//skip a comment line used as header in parameterfile
			runlog << dummystring << " \n";
			int SavePropertyLog; //saving property log file, useful if calibrating any of the skytype scores (ok for long runs)
				parameterfile >> SavePropertyLog; getline(parameterfile,dummystring);
				runlog << "\t" << SavePropertyLog << dummystring << " \n"; 
			double HalfWidth;  // half width of broadening Gaussians
				parameterfile >> HalfWidth; getline(parameterfile,dummystring);
				runlog << "\t" << HalfWidth << dummystring << " \n"; 
		//SkyTypePrms, they should be copied/pasted from the appropriate sheet in Haloloop-Skytype-HaloScore-Parametercollection.xlsx
			double SkyTypePrmsM[4][10]; // the vector for average values for each sky type
			//	index		meaning
			//	1	[4]		sky type index: 0-CS, 1-PCL, 2-CLD, 3-CLR
			//	2	[10]	property index:
			//					0-	B	 slope
			//					1-	B	 intercept
			//					2-	B	 ARIvar
			//					3-	B	 CV
			//					4-	G	 slope
			//					5-	G	 intercept
			//					6-	G	 ARIvar
			//					7-	R	 slope
			//					8-	R	 intercept
			//					9-	R	 ARIvar
			double SkyTypePrmsICM [4][10][10]; //Inverted Covariance Matrix for each sky type
			// index interpretation analogous to above.
			getline(parameterfile,dummystring);// take in the header line for skytype table
			runlog << dummystring << " \n"; // write it to runlog
			getline(parameterfile,dummystring);//take the line containing the generation date of the skytype prms
			runlog << dummystring << " \n";//write it to runlog
			for (int itype = 0; itype <4; itype++)
			{
				getline(parameterfile,dummystring);// the number of records used for the current skytype averages
				runlog << dummystring << " \n"; // write it to runlog
				getline(parameterfile,dummystring);// M vector line
				runlog << dummystring << " \n"; // write it to runlog
				for (int iprop = 0; iprop < 10; iprop++ )
				{ 
					parameterfile >> SkyTypePrmsM[itype][iprop];//reading the M vector
					runlog << SkyTypePrmsM[itype][iprop] << "\t";
				} 
				runlog << "\n";
				getline(parameterfile,dummystring);// this is here to get to eol after matrix was read
				getline(parameterfile,dummystring);// inverse matrix line
				runlog << dummystring << " \n"; // write it to runlog
				for (int iprop = 0; iprop < 10; iprop++ )
				{ 
					for (int jprop = 0; jprop < 10; jprop++)
					{
						parameterfile >> SkyTypePrmsICM [itype][iprop][jprop]; // reading the inverse covariance matrix
						runlog << SkyTypePrmsICM [itype][iprop][jprop] << "\t" ; 
					}//done jprop
					runlog << "\n";
				}// done iprop
				getline(parameterfile,dummystring);// this is here to get to eol after matrix was read
			}//done itype
		//HaloScorePrms, they should be copied/pasted from the appropriate sheet in Haloloop-Skytype-HaloScore-Parametercollection.xlsx
			double HaloScorePrmsM[31]; // the vector for average values for the halo score
			//	index		meaning
			//	1	[31]	property index:
			//				0	B slope				1	B intercept			2	B ARIvar			3	CV				4	B LocMaxDer
			//				5	B LocPMDer			6	B LocMinDer			7	B ValMaxDer			8	B ValMinDer		9	B nmax
			//				10	B stdLMaD			11	B stdLPMD			12	B stdLMiD			13	G slope			14	G intercept
			//				15	G ARIvar			16	G LocMaxDer			17	G LocPMDer			18	G LocMinDer		19	G ValMaxDer
			//				20	G ValMinDer			21	G nmax				22	R slope				23	R intercept		24	R ARIvar
			//				25	R LocMaxDer			26	R LocPMDer			27	R LocMinDer			28	R ValMaxDer		29	R ValMinDer
			//				30	R nmax							

			double HaloScorePrmsICM [31][31]; //Inverted Covariance Matrix for halo score
			// index interpretation analogous to above.
			getline(parameterfile,dummystring);// take in the header line for haloscore table
			runlog << dummystring << " \n"; // write it to runlog
			getline(parameterfile,dummystring);//take the line containing the generation date of the haloscore prms
			runlog << dummystring << " \n";//write it to runlog
			getline(parameterfile,dummystring);// the number of records used for the current haloscore averages
			runlog << dummystring << " \n"; // write it to runlog
			getline(parameterfile,dummystring);// M vector line
			runlog << dummystring << " \n"; // write it to runlog
			for (int iprop = 0; iprop < 31; iprop++ )
			{ 
				parameterfile >> HaloScorePrmsM[iprop];//reading the M vector for halo score
				runlog << HaloScorePrmsM[iprop] << "\t";
			} 
			runlog << "\n";
			getline(parameterfile,dummystring);// this is here to get to eol after matrix was read
			getline(parameterfile,dummystring);// inverse matrix line
			runlog << dummystring << " \n"; // write it to runlog
			for (int iprop = 0; iprop < 31; iprop++ )
			{ 
				for (int jprop = 0; jprop < 31; jprop++)
				{
					parameterfile >> HaloScorePrmsICM[iprop][jprop]; // reading the inverse covariance matrix
					runlog <<  HaloScorePrmsICM[iprop][jprop] << "\t"; // reading the inverse covariance matrix
				}//done jprop
				runlog << "\n";
			}// done iprop	
			getline(parameterfile,dummystring);// this is here to get to eol after matrix was read
			
		parameterfile.close();
		runlog << "** Finished reading "<< parameterfilename << "\n";
	//*****************done Reading the input parameter file **********************************************
				
	//***************  writing some output file headers ************************
		string linefitfilename;
		sc << imagelistfile << "-haloloop-out-LineFitPrms.csv"; linefitfilename = sc.str(); sc.str("");//assign the name of the parameterfile
		ofstream linefitfile;
		if (SaveLineFitPRMs == 1)
		{
			linefitfile.open(linefitfilename);
			runlog << "** saving the line fit parameters for each image in " << linefitfilename << " \n";
			linefitfile << "rfitlimmin," << rfitlimmin  << ", rfitlimmax," << rfitlimmax << endl;
			linefitfile << "image, TR-B-slope, TR-B-intercept, TR-B-R2, TR-G-slope, TR-G-intercept, TR-G-R2, TR-R-slope, TR-R-intercept, TR-R-R2,  ";
			linefitfile << "BR-B-slope, BR-B-intercept, BR-B-R2, BR-G-slope, BR-G-intercept, BR-G-R2, BR-R-slope, BR-R-intercept, BR-R-R2,   ";
			linefitfile << "BL-B-slope, BL-B-intercept, BL-B-R2, BL-G-slope, BL-G-intercept, BL-G-R2, BL-R-slope, BL-R-intercept, BL-R-R2,  ";
			linefitfile << "TL-B-slope, TL-B-intercept, TL-B-R2, TL-G-slope, TL-G-intercept, TL-G-R2, TL-R-slope, TL-R-intercept, TL-R-R2 \n";
		}
		
		string SkyTypeFilename;
		sc << imagelistfile << "-haloloop-out-SkyTypeLog.csv"; SkyTypeFilename = sc.str(); sc.str("");//assign the name of the parameterfile
		ofstream SkyTypeFile;
		SkyTypeFile.open(SkyTypeFilename);
		runlog << "** saving skytype probabilities to " << SkyTypeFilename << " \n";
		SkyTypeFile << " origin, UTZ, image,AzimuthAngle,ZenithAngle,, ";
		SkyTypeFile << "All-DomSky, Halo?,,  TR-HaloDos, BR-HaloDos, BL-HaloDos, TL-HaloDos,All-HaloDos,, All-CS, All-PCL, All-CLD, All-CLR , All-N/A, All-HaloScore,,";
		SkyTypeFile << "TR-DomSky, TR-CS, TR-PCL, TR-CLD, TR-CLR ,TR-N/A,TR-HaloScore,TR-D,,";
		SkyTypeFile << "BR-DomSky, BR-CS, BR-PCL, BR-CLD, BR-CLR ,BR-N/A,BR-HaloScore,BR-D,,";
		SkyTypeFile << "BL-DomSky, BL-CS, BL-PCL, BL-CLD, BL-CLR ,BL-N/A,BL-HaloScore,BL-D,,";
		SkyTypeFile << "TL-DomSky, TL-CS, TL-PCL, TL-CLD, TL-CLR ,TL-N/A,TL-HaloScore,TL-D  \n";
		
		string HaloLogFilename;
		sc << imagelistfile << "-haloloop-out-HaloLog.csv"; HaloLogFilename = sc.str(); sc.str("");//assign the name of the parameterfile
		ofstream HaloLogFile;
		HaloLogFile.open(HaloLogFilename);
		runlog << "** saving Halo notes to " << HaloLogFilename << " \n";
		HaloLogFile << " origin, UTZ, image,AzimuthAngle,ZenithAngle,, ";
		HaloLogFile << "All-DomSky, Halo?, All-CS, All-PCL, All-CLD, All-CLR, All-N/A, All-HaloScore,,";
		HaloLogFile << "TR-DomSky, TR-CS, TR-PCL, TR-CLD, TR-CLR ,TR-N/A, TR-HaloScore,TR-D,,";
		HaloLogFile << "BR-DomSky, BR-CS, BR-PCL, BR-CLD, BR-CLR ,BR-N/A,BR-HaloScore,BR-D,,";
		HaloLogFile << "BL-DomSky, BL-CS, BL-PCL, BL-CLD, BL-CLR ,BL-N/A,BL-HaloScore,BL-D,,";
		HaloLogFile << "TL-DomSky, TL-CS, TL-PCL, TL-CLD, TL-CLR ,TL-N/A,TL-HaloScore,TL-D \n";


		string PropertyLogFile;
		sc << imagelistfile << "-haloloop-out-PropertyLog.csv"; PropertyLogFile = sc.str(); sc.str("");//assign the name of the parameterfile
		ofstream PropertyLog;
		if (SavePropertyLog == 1)
		{
			PropertyLog.open(PropertyLogFile);
			runlog << "** saving scored properties to  " << PropertyLogFile << " ******************************** \n";
			PropertyLog << "******************  " << PropertyLogFile << "  ************************************  \n ";
			PropertyLog << "****lines from this table can be directly copied and pasted to the parameter collection file **  \n ";
			PropertyLog << "****              haloloop/src/Halolopp-SkyType-HaloScore-Parametercollection.xlsx       ****  \n \n";
            PropertyLog << " image, quadrant, B,,,,,,,,,,,,,G,,,,,,,,,,,,,R,,,,,,,,,,,,,B,B,G,G,R,R,B,G,R \n ,,";			
			for (int i = 0; i <3; i++) { PropertyLog << " slope, intercept, ARIvar, CV, LocMaxDer, LocPMDer, LocMinDer, ValMaxDer, ValMinDer, nmax, stdLMaD, stdLPMD, stdLMiD, ";}
			PropertyLog << " mxtozero,maxtomin,mxtozero,maxtomin,mxtozero,maxtomin,DerRange,DerRange,DerRange, \n";			
		}	 
	// ************ done writing writing some output file headers ****
	
	// initiate a time keeper variable, which triggers an empty line in skytype log if a time gap in the data is longer than 6 hours (nights, malfunctions)
	double PreviousTime = 0.0;
	int ContImage = 0; // a counter that keeps increasing as long as the image time is continously increasing, and restarts if a gap occurs
	
	//initiate the buffer vectors, which will keep preceding and following skytype out put in order to allow a dos construction
	string BufferString[200][2];
	// 		[200]	ibuff		the number of records in the buffer with a range of 3w, containing 6w+1 records. Max width for dos is then 33 (should never be necessary)
	//							the current record shold always sit at ibuff = 3w
	//		[2]		istr		0 = the skytypelog output string before the dos, 1 = the skytypelog output string after the dos
	double BufferDos[200][5];	
	// 		[200]	ibuff		the number of records in the buffer with a range of 3w, containing 6w+1 records. Max width for dos is then 33 (should never be necessary)
	//							the current record shold always sit at ibuff = 3w
	//		[5]		idos		0 = TR
	//							1 = BR
	//							2 = BL
	//							3 = TL
	//							4 = all sky		
	// Define the buffer properties
	int BuffRange = 3 * int(HalfWidth+0.5); // to include up to the 3sigma range of a bell curve
	int BuffTop = 100 + BuffRange;
	int BuffBot = 100 - BuffRange;

	
	
	//******************************Start processing images *******************************************************
		runlog << "******************** Starting Images ***********************************\n";
		// open the temporary processlist as an input stream
		ifstream processlistin("haloloop-processlist.txt");
		//begin the loop over all images
		for (imagecounter = 0; imagecounter < totalimagecount; imagecounter++) 
		{
			//import the information for the current image from the temporary process list
			string origin;
			int year,month,day,hour,minute,second;
			processlistin >> CurrentImage >> origin >> year >> month >> day >> hour >> minute >> second;
			
			
			//Load the current image to SourceImageMat
			Mat SourceImageMat = imread(CurrentImage, 1); 
			
			//a check that the CurrentImage is existent
			if (!SourceImageMat.data)
			{ runlog << "Could not find " << CurrentImage << endl;} //if the image is empty no procedures are applied
			else {//executed if image not empty 
				//**** image dimensions as used throughout the process for the current image
				int cols = SourceImageMat.cols, rows = SourceImageMat.rows;
				
				if (imagecounter%500 == 0) {cout   << " \n "<< imagecounter ;}
				runlog << " working on " << CurrentImage  << "\n";

				// ****************apply color corrections *******************************************
				//find the average color for a small rectangle above center
				Mat SourceColorMat; //draw the ellipse as specified by the image pre-processing parameters
				SourceImageMat.copyTo(SourceColorMat);//initialize with the original image
				int avcol[3]; avcol[0]=0; avcol[1]=0; avcol[2]=0;
				int counter = 0;
				int startcol = int(0.6 * double(cols)), endcol = int(0.8 * double(cols));
				int startrow = int(0.3* double(rows)), endrow = int(0.45 * double(rows));
				for (int icol = startcol; icol < endcol; icol++)
				{
					for (int irow = startrow; irow < endrow; irow++)
					{
						Vec3b bgr = SourceImageMat.at<Vec3b>(irow,icol);
						avcol[0] += bgr[0];
						avcol[1] += bgr[1];
						avcol[2] += bgr[2];
						counter += 1;
					}
				}
				avcol[0] = avcol[0] /counter;avcol[1] = avcol[1] /counter;avcol[2] = avcol[2] /counter;
				//cout << "average color: " << avcol[0] << "\t" << avcol[1] << "\t" << avcol[2] << "\n";
				double TargetColor[3];
				TargetColor[0] = RelcolShift[0]*double(avcol[0]);if (TargetColor[0] > 255.0) TargetColor[0] = 255.0;
				TargetColor[1] = RelcolShift[1]*double(avcol[1]);if (TargetColor[1] > 255.0) TargetColor[1] = 255.0;
				TargetColor[2] = RelcolShift[2]*double(avcol[2]);if (TargetColor[2] > 255.0) TargetColor[2] = 255.0;
				double alpha = RelcolShift[3];
				for(int icol = 0; icol < cols; icol++)
				{
					for (int irow = 0; irow < rows; irow++ )
					{
						Vec3b CC = SourceImageMat.at<Vec3b>(irow, icol);
						double B,G,R; B = double(CC[0]); G = double(CC[1]); R = double(CC[2]);
						double ccsum1 = B + G + R;
						B = B + (TargetColor[0] - B) * alpha; 
						G = G + (TargetColor[1] - G) * alpha; 
						R = R + (TargetColor[2] - R) * alpha; 
						double ccsum2 = B + G + R;
						B = B * ccsum1 / ccsum2; if (B > 255.0) B = 255.0; if (B < 0) B = 0.0;
						G = G * ccsum1 / ccsum2; if (G > 255.0) G = 255.0; if (G < 0) G = 0.0;
						R = R * ccsum1 / ccsum2; if (R > 255.0) R = 255.0; if (R < 0) R = 0.0;
						Vec3b FC; FC[0] = int(B); FC[1] = int(G); FC[2] = int(R);
						SourceColorMat.at<Vec3b>(irow, icol)= FC;
					}// done irow
				}// done icol									
				SourceColorMat.copyTo(SourceImageMat);//initialize with the original image
				// ***************done apply color corrections *******************************************	

				
				//*** identify the radius and center position of the horizon ellipse for the image in the pixel system
				Point2f HorizonCenter; HorizonCenter.x = double(cols)/2.0+HCenterShift.x; HorizonCenter.y = double(rows)/2.0+HCenterShift.y; //contains x and y coordinates of the center of the ellipse, these are initial values
				Point2f HorizonAxes; HorizonAxes.x=double(cols)/2.0*HXadjust; HorizonAxes.y=HorizonAxes.x*eps; //contains the semiaxes of the horizon ellipse, initial value vaLues for both axes
				if (SaveHorizonMat == 1) //this is an option set in the header, avoid with big numbers of files
				{
					Mat SourceHorizonMat; //draw the ellipse as specified by the image pre-processing parameters
					SourceImageMat.copyTo(SourceHorizonMat);//initialize with the original image
					ellipse( SourceHorizonMat, HorizonCenter, Size( HorizonAxes.x, HorizonAxes.y), 0, 0, 360, Scalar( 0, 0, 255), 1, 8, 0);
					circle( SourceHorizonMat, HorizonCenter, 6, Scalar( 0, 0, 255), -1, 8, 0); //the center
					stringstream sc;sc.str(""); sc << CurrentImage<<"-1-HorizonMat.jpg"; string filename = sc.str(); imwrite(filename, SourceHorizonMat);
					runlog << "\t saving source with marked horizon in " << filename << "\n";
				}
								

			
				//********Removing the distortion ********************************************
				// stretch the ellipse to a circular view and move from zenith projection to arc projection, all in pixel coordinates (see documentation for details)
				double HorizonRadius = HorizonAxes.x; //assigning the radius of the circular horizon port
				Mat CircImageMat; //Mat to only contain the source image stretched into a circle
				SourceImageMat.copyTo(CircImageMat);//will contain the circular image, stretched from the elliptical view
				{  // begin transformation from elliptic horizon to circular horizon
					double vcx = HorizonCenter.x, vcy = HorizonCenter.y; //grab center
					for (int circrow = 0; circrow < rows; circrow++){
						for (int circcol = 0; circcol < cols; circcol++){
							Vec3b srcbgr; //a vector to hold the intensity values in the source image
							int sourcecol, sourcerow; //pixel coord in source image
							//find the jrow and jcol in the source image that belongs to irow and icol in CircImageMat
							double xcirc, ycirc, Rxcirc, Rycirc, R;
							xcirc = double(circcol); //column coordinate in circimg
							ycirc = double(circrow); //row coordinate in circimg
							Rxcirc = xcirc - vcx; Rycirc = ycirc - vcy; //coordinate relative to horizon center
							double R2 = Rxcirc*Rxcirc + Rycirc*Rycirc; R = sqrt(R2);// R is the vector from center of view to current pixel in circ image							
							double Rxsource = R2 / (eps*eps*Rxcirc*Rxcirc + Rycirc*Rycirc); Rxsource = eps*Rxcirc*sqrt(Rxsource);//ellipse transformation, described in documentation
							double Rysource = Rycirc*Rxsource / Rxcirc;
							if (Rycirc == 0) { Rxsource = Rxcirc; Rysource = 0; }
							if (Rxcirc == 0){ Rxsource = 0; Rysource = eps*Rycirc; }
							double ysource = vcy + Rysource, xsource = vcx + Rxsource; // position of the pixel in the SourceImageMat
							sourcecol = int(xsource + 0.5); sourcerow = int(ysource + 0.5); //assigned integer pixel position in zenith view				
							if (sourcecol < 0) {sourcecol = 0;}
							if (sourcerow < 0) {sourcerow = 0;} //two lines to deal with edge pixels
							if (sourcecol >= cols) {sourcecol = cols - 1;} 
							if (sourcerow >= rows) {sourcerow = rows - 1;}
							if ( R >= HorizonRadius-0.5 )
							{
								srcbgr[0]=0;srcbgr[1]=0;srcbgr[2]=0;//blackening the area outside of view circle
							}
							else
							{
								srcbgr = SourceImageMat.at<Vec3b>(sourcerow, sourcecol);//retrieve the bgr values from source image
							}
							CircImageMat.at<Vec3b>(circrow, circcol)= srcbgr;//assign color of pixel in circular image
						}//circcol
					}//circrow
					if (SaveCircImageMat == 1) //this is an option set in the header, avoid with big numbers of files
					{
						Mat CircTempMat; //internal name for the circ image with horizon circle drawn
						CircImageMat.copyTo(CircTempMat);//initialize with the original image
						circle( CircTempMat, HorizonCenter, HorizonRadius, Scalar( 0, 0, 255), 1, 8, 0); //the horizon circle
						circle( CircTempMat, HorizonCenter, 6, Scalar( 0, 0, 255), -1, 8, 0); //the center
						stringstream sc;sc.str(""); sc << CurrentImage<<"-2-CircImage.jpg"; string filename = sc.str(); imwrite(filename, CircTempMat);
						runlog << "\t saving CircImageMat in " << filename << "\n";
					}//end SaveCircMat
				}// end transformation from elliptic horizon to circular horizon

				
				//*********************** correct the north alignment
				Mat NorthAlignMat; //Mat to hold the rotation-calibrated image, worked in pixel coordinates
				CircImageMat.copyTo(NorthAlignMat);
				if (NorthAlignAngle != 0.0)
				{
					double angle = NorthAlignAngle * Pi / 180.0; //angle in radians
					double vcx = HorizonCenter.x, vcy = HorizonCenter.y; //grab center
					for (int narow = 0; narow < rows; narow++){
						for (int nacol = 0; nacol < cols; nacol++){
							Vec3b circbgr; //a vector to hold the intensity values in the source image
							int circcol, circrow; //pixel coord in source image
							//find the circrow and circcol in the CircImageMat that belongs to irow and icol in NorthAlignMat
							double xna, yna, Rxna, Ryna;
							xna = double(nacol); //column coordinate in NorthAlignMat
							yna = double(narow); //row coordinate in NorthAlignMat
							Rxna = xna - vcx; Ryna = yna - vcy; //coordinate relative to horizon center	
							double R = Rxna * Rxna + Ryna * Ryna; R = sqrt(R);
							if (R < HorizonRadius)
							{
								double Rxcirc =  Rxna*cos(angle) + Ryna*sin(angle);
								double Rycirc = -Rxna*sin(angle) + Ryna*cos(angle);
								double ycirc = vcy + Rycirc, xcirc = vcx + Rxcirc; // position of the pixel in the CircImageMat
								circcol = int(xcirc + 0.5); circrow = int(ycirc + 0.5); //assigned integer pixel position in zenith view				
								if (circcol < 0) {circcol = 0;}
								if (circrow < 0) {circrow = 0;} //two lines to deal with edge pixels
								if (circcol >= cols) {circcol = cols - 1;} 
								if (circrow >= rows) {circrow = rows - 1;}	
							} else
							{
								circrow = narow;
								circcol = nacol;
							}
							circbgr = CircImageMat.at<Vec3b>(circrow, circcol);//retrieve the bgr values from source image
							NorthAlignMat.at<Vec3b>(narow, nacol)= circbgr;//assign color of pixel in circular image
						}//nacol
					}//narow
					if (SaveNorthAlignMat == 1) //this is an option set in the header, avoid with big numbers of files
					{
						Mat NATempMat; //internal name for the circ image with horizon circle drawn
						NorthAlignMat.copyTo(NATempMat);//initialize with the original image
						circle( NATempMat, HorizonCenter, HorizonRadius, Scalar( 0, 0, 255), 3, 8, 0); //the horizon circle
						circle( NATempMat, HorizonCenter, 6, Scalar( 0, 0, 255), -1, 8, 0); //the center
						stringstream sc;sc.str(""); sc << CurrentImage <<"-3-NAligned.jpg"; string filename = sc.str(); 
						imwrite(filename, NATempMat);
						runlog << "\t saving NorthAlignMat in " << filename << "\n";
					}//end SaveNorthAlignMat
				}//*********************** done correct the north alignment
				
												
				//**** Solar position in image ****************************************************************
				Point2f SunCenter; //in pixel coordinates
				double AzimuthAngle, ZenithAngle;
				{ //getting zenith angle and azimuth
					double Latitude, Longitude; //coordinates of ARM TSI site
					if (origin == "nsa" )
					{  // 71° 19′ 22.8″ N, 156° 36′ 32.4″ W					
						Latitude = 71.0 + 19.0/60.0 + 22.8/3600; //local latitude in degrees
						Longitude = -156.0 - 36.0/60.0 - 32.4/3600; //local longitude in degrees
					}				
					if (origin == "sgp" )
					{  // 36° 36′ 18″ N, 97° 29′ 6″ W
						Latitude = 36.0 + 36.0/60.0 + 18.0/3600; //local latitude in degrees
						Longitude = -97.0 - 29.0/60.0 - 6.0/3600; //local longitude in degrees
					}
					if (origin == "ena" )
					{  // 39° 5′ 29.76″ N, 28° 1′ 32.52″ W
						Latitude = 39.0 + 5.0/60.0 + 29.76/3600; //local latitude in degrees
						Longitude = -28.0 - 1.0/60.0 - 32.52/3600; //local longitude in degrees
					}
					if (origin == "oli" )
					{  // 70deg  29’ 42’’ N, 149deg 53’ 9.6’’ W
						Latitude = 70.0 + 29.0/60.0 + 42.0/3600; //local latitude in degrees
						Longitude = -149.0 - 53.0/60.0 - 9.6/3600; //local longitude in degrees
					}					
					//get the date as the day of the year, found at https://stackoverflow.com/questions/19110675/calculating-day-of-year-from-date
					tm timeinfo = {}; timeinfo.tm_year = year - 1900;
					timeinfo.tm_mon = month - 1; timeinfo.tm_mday = day; mktime(&timeinfo);
					double YearDay = timeinfo.tm_yday; //whole number of days at the beginning of current day
					double YearLength = 365.25; // length of year
					double VernalEquinox = 79.0; //day of March 20, 5pm, 2019
					double Obliquity = 23.4;// in units of degrees
					//times
					double GSTdeg = ( YearDay - VernalEquinox )/YearLength*360.0 + 180.0; if (GSTdeg > 360.0 ){GSTdeg = GSTdeg -360.0;} //Greenwich Sidereal Time in degrees
					double LSTdeg = GSTdeg + Longitude + (double(hour)+double(minute)/60.0 + double(second)/3600.0)*360.0/24.0; //local sideric time
					if (LSTdeg > 360.0){LSTdeg = LSTdeg -360.0;}
					//ephemerides
					double Declination = Obliquity * sin(( YearDay + double(hour)/24.0 + double(minute)/1440.0 - VernalEquinox)/YearLength*2.0*Pi); //declination on this day
					double RightAscension = GSTdeg +180.0; if (RightAscension > 360.0){RightAscension = RightAscension -360.0;}//Ascension in degrees
					//horizon coordinates
					double HourAngle = LSTdeg - RightAscension; if (HourAngle > 360.0){HourAngle = HourAngle -360.0;}
					//AzimuthAngle					
					double HourAnglerad = HourAngle * Pi/180.0; 
					double SinHour = sin( HourAnglerad ), CosHour = cos(HourAnglerad);
					double Latituderad = Latitude * Pi/180.0, SinLat = sin(Latituderad), CosLat = cos(Latituderad);
					double Declinationrad = Declination * Pi/180.0, CosDec = cos(Declinationrad), SinDec = sin(Declinationrad);
					//azimuth angle using the transformation given in https://en.wikipedia.org/wiki/Celestial_coordinate_system
					double xTanArg = -SinLat * CosDec * CosHour + CosLat * SinDec;
					double yTanArg = CosDec * SinHour;
					AzimuthAngle = - atan2(yTanArg, xTanArg )*180.0/Pi;					
					//elevation
					double Sina = SinLat*SinDec + CosLat * CosDec * CosHour;
					ZenithAngle = 90.0 - asin(Sina)*180.0/Pi;
					SunCenter.x = HorizonCenter.x + ZenithAngle / MaximumZenithAngle * HorizonRadius * sin( AzimuthAngle*Pi/180.0);
					SunCenter.y = HorizonCenter.y - ZenithAngle / MaximumZenithAngle * HorizonRadius * cos( AzimuthAngle*Pi/180.0);	
					if (SaveSunMarked == 1) //this is an option set in the header, avoid with big numbers of files
					{
						Mat SunTempMat; //internal name for the circ image with horizon circle drawn
						//***** entered 6-25-2018: I do not need the distortion removal! Use the circimage directly.
						NorthAlignMat.copyTo(SunTempMat);//initialize with the original image
						circle( SunTempMat, HorizonCenter, HorizonRadius, Scalar( 0, 0, 255), 1, 8, 0); //the horizon circle
						circle( SunTempMat, HorizonCenter, 6, Scalar( 0, 0, 255), -1, 8, 0); //the center
						circle( SunTempMat, SunCenter, 6, Scalar( 0, 255, 255), -1, 10, 0); //the sun
						stringstream sc;sc.str(""); sc << CurrentImage<<"-4-SunMarked.jpg"; string filename = sc.str(); imwrite(filename, SunTempMat);
						runlog << "\t saving CircImageMat with sun in " << filename << "\n";
					}//end SaveArcMat	
					//runlog << "\t \t AzimuthAngle (deg): \t" << AzimuthAngle << "\t \t ZenithAngle(deg): \t" << ZenithAngle << endl;
				} // done zenithangle and azimuth								
				//**** Done Solar position in image ***********************************************************

				//******** Removing the perspective distortion ********************************************
				Mat ArcImageMat;//this is the arc projection of the masked circular image, all worked in pixel coordinates
				NorthAlignMat.copyTo(ArcImageMat);//this will contain the arc projection of the circular image
				{//begin zenith projection to arc projection
					double thetaH = MaximumZenithAngle*Pi/180.0; //zenith angle of the horizon circle
					double sinthH = sin(thetaH);
					double hcx = HorizonCenter.x, hcy = HorizonCenter.y;
					for (int arcrow = 0; arcrow < rows; arcrow++){
						for (int arccol = 0; arccol < cols; arccol++){
							int circcol, circrow; //pixel coord in zenith image
							//find the circcol and cricrow in the zenith image (circimagemat) that belongs to arccol and arcrow in the new arc image
							double xarc, yarc, sx, sy, s;
							xarc = double(arccol); //column coordinate in arcimg
							yarc = double(arcrow); //row coordinate in arcimg
							sx = xarc - hcx;//vector from center to current pixel
							sy = yarc - hcy;
							s = sx*sx + sy*sy; 
							s = sqrt(s);// s is the vector from center of view to current pixel in arc image	
							if (s < HorizonRadius)
							{
								if (abs(sx) < 0.5) {sx = 0.0;} 
								if (abs(sy) < 0.5) {sy = 0.0;} //avoid overflow if components are too small
								double sex = sx / s, sey = sy / s;//unit vector for s
								double r = HorizonRadius * sin( s * thetaH / HorizonRadius);
								r = r / sinthH;//length of a in zenith img		
								if (r > HorizonRadius) {r = HorizonRadius;} //avoid bleading the mask into the image
								double xcirc = hcx + r*sex; // position of the pixel in the zenith image
								double ycirc = hcy + r*sey;
								circcol = int(xcirc + 0.5);
								circrow = int(ycirc + 0.5); //assigned integer pixel position in zenith view				
								if (circcol < 0) { circcol = 0;}
								if (circrow < 0) { circrow = 0;} //two lines to deal with edge pixels
								if (circcol >= cols) {circcol = cols - 1;}
								if (circrow >= rows) {circrow = rows - 1;}	
							}
							else
							{
								circrow = arcrow; circcol = arccol;
							}							
							ArcImageMat.at<Vec3b>(arcrow, arccol)=NorthAlignMat.at<Vec3b>(circrow, circcol);
						}//arccol
					}//arcrow
					//adjusting solar position
					Point2f relsunpos = SunCenter - HorizonCenter; double r = norm(relsunpos);
					//double r = relsunpos.x*relsunpos.x+relsunpos.y*relsunpos.y; r = sqrt(r);
					double ss=r/HorizonRadius*sinthH;
					if (ss > 1.0 ) {ss = 1.0;}
					double Suns = HorizonRadius/thetaH * asin(ss);
					if (Suns > HorizonRadius) {Suns = HorizonRadius-0.5;} //for very low sun situations, the placement can be outside of horizon cirlce. This avoids that.
					SunCenter.x = HorizonCenter.x + Suns * sin( AzimuthAngle*Pi/180.0);
					SunCenter.y = HorizonCenter.y - Suns * cos( AzimuthAngle*Pi/180.0);					
					if (SaveArcImageMat == 1) //this is an option set in the header, avoid with big numbers of files
					{
						Mat CircTempMat; //internal name for the circ image with horizon circle drawn
						ArcImageMat.copyTo(CircTempMat);//initialize with the original image
						circle( CircTempMat, HorizonCenter, HorizonRadius, Scalar( 0, 0, 255), 1, 8, 0); //the horizon circle
						circle( CircTempMat, HorizonCenter, 6, Scalar( 0, 0, 255), -1, 8, 0); //the center
						circle( CircTempMat, SunCenter, 6, Scalar( 0, 255, 255), -1, 8, 0); //the sun
						stringstream sc;sc.str(""); sc << CurrentImage<<"-5-ArcImage.jpg"; string filename = sc.str(); imwrite(filename, CircTempMat);
						runlog << "\t saving ArcImageMat in " << filename << "\n";
					}//end SaveCircMat
				}//done zenith projection to arc projection
				//******** Done Removing the perspective distortion ********************************************				

				
				//***** mask the shadow strip and the cameramount before starting halo analysis **********************
				Mat SkyMaskMat; //internal name for image with features other than sky masked, also worked in pixels
				ArcImageMat.copyTo(SkyMaskMat); // this step can not be skipped; it assigns the size of SkyMaskMat which allows the color assignments below
				{ // begin shadow mask
					Point2f VecCentSun = SunCenter - HorizonCenter; VecCentSun = VecCentSun / norm(VecCentSun); //radial vector center through sun
					//camera mount https://stackoverflow.com/questions/40284369/opencv-fillconvexpoly-function-in-c-throws-exception, order of points matters
					vector<Point> MountPoints;
					MountPoints.push_back(Point(int(HorizonCenter.x-1), int(HorizonCenter.y+2))); //left lower corner
					MountPoints.push_back(Point(int(HorizonCenter.x+1), int(HorizonCenter.y+2))); //right lower corner
					MountPoints.push_back(Point(int(HorizonCenter.x+4), int(HorizonCenter.y-HorizonRadius-5))); //right upper corner
					MountPoints.push_back(Point(int(HorizonCenter.x-4), int(HorizonCenter.y-HorizonRadius-5))); //left upper corner
					fillConvexPoly(SkyMaskMat,MountPoints,Scalar(0,0,0));
					//shadow strip
					Point2f etan; etan.x = VecCentSun.y; etan.y = -VecCentSun.x; //tangential vector, clockwise
					vector<Point> ShadePoints;
					int ixa = int(HorizonCenter.x + (0.5*ShadowStripWidth+ShadowStripOffset)*etan.x + 0.5);
					int iya = int(HorizonCenter.y + (0.5*ShadowStripWidth+ShadowStripOffset)*etan.y + 0.5);
					int ixb = int(HorizonCenter.x + 1.1*HorizonRadius*VecCentSun.x + (0.5*ShadowStripWidth+ShadowStripOffset)*etan.x + 0.5);
					int iyb = int(HorizonCenter.y + 1.1*HorizonRadius*VecCentSun.y + (0.5*ShadowStripWidth+ShadowStripOffset)*etan.y + 0.5);
					int ixc = int(HorizonCenter.x + 1.1*HorizonRadius*VecCentSun.x - (0.5*ShadowStripWidth-ShadowStripOffset)*etan.x + 0.5);
					int iyc = int(HorizonCenter.y + 1.1*HorizonRadius*VecCentSun.y - (0.5*ShadowStripWidth-ShadowStripOffset)*etan.y + 0.5);
					int ixd = int(HorizonCenter.x - (0.5*ShadowStripWidth-ShadowStripOffset)*etan.x + 0.5);
					int iyd = int(HorizonCenter.y - (0.5*ShadowStripWidth-ShadowStripOffset)*etan.y + 0.5);
					ShadePoints.push_back(Point(ixa, iya));
					ShadePoints.push_back(Point(ixb, iyb));
					ShadePoints.push_back(Point(ixc, iyc));
					ShadePoints.push_back(Point(ixd, iyd));
					fillConvexPoly(SkyMaskMat,ShadePoints,Scalar(0,0,0));	
					circle( SkyMaskMat, HorizonCenter, ShadowStripWidth/2.0, Scalar( 0, 0, 0), -1, 10, 0); //the sun					
					if (SaveSkyMaskMat == 1) //this is an option set in the header, avoid with big numbers of files
					{
						Mat SunTempMat; 
						SkyMaskMat.copyTo(SunTempMat);//initialize with the original image
						circle( SunTempMat, HorizonCenter, HorizonRadius, Scalar( 0, 0, 255), 1, 8, 0); //the horizon circle
						circle( SunTempMat, HorizonCenter, 6, Scalar( 0, 0, 255), -1, 8, 0); //the center
						circle( SunTempMat, SunCenter, 6, Scalar( 0, 255, 255), -1, 10, 0); //the sun
						stringstream sc;sc.str(""); sc << CurrentImage<<"-6-SkyMaskMat.jpg"; string filename = sc.str(); imwrite(filename, SunTempMat);
						runlog << "\t saving SkyMaskMat in " << filename << "\n";
					}//end SaveSkyMaskMat					
				} //done shadow mask
				//***** done mask the shadow strip and the cameramount before starting halo analysis **********************
		
				
				//****************** create the Local Sky Map (SRM) ******************************************************
				//**** that is a rotated (SkyRotMat) and subsequently cropped to 150 by 150 around sun (LocalSkyMat) version of the masked image
				double LSMCropSize = HorizonRadius / MaximumZenithAngle * LSMCropSizesd; //ascertaining crop size in pixels
				Mat LocalSkyMat = Mat::zeros(LSMCropSize,LSMCropSize, CV_8UC3);
				{
					Mat SkyRotMat; //will contain the rotated version of tne masked sky image
					SkyMaskMat.copyTo(SkyRotMat);
					double gamma = AzimuthAngle-180.0; gamma = gamma * Pi/180.0; //rotation angles of sun into south
					double cosgamma = cos(gamma), singamma = sin(gamma);
					double hcx = HorizonCenter.x, hcy = HorizonCenter.y, scx = SunCenter.x, scy = SunCenter.y;
					for (int SRMrow = 0; SRMrow < rows; SRMrow++){
						for (int SRMcol = 0; SRMcol < cols; SRMcol++){
							Vec3b SRMbgr;
							double xSMM = hcx + (double(SRMcol)-hcx)*cosgamma - (double(SRMrow)-hcy)*singamma;
							double ySMM = hcy + (double(SRMcol)-hcx)*singamma + (double(SRMrow)-hcy)*cosgamma;
							int SMMrow = int( ySMM + 0.5); int SMMcol = int( xSMM +0.5 );
							if ((SMMrow >= 0) && (SMMrow < rows) && (SMMcol >= 0) && (SMMcol < cols))
							{
								SRMbgr= SkyMaskMat.at<Vec3b>(SMMrow, SMMcol);//retrieve the bgr values from SMM												
							} else 
							{
								SRMbgr[0]=0; SRMbgr[1]=0, SRMbgr[2]=0;
							}
							SkyRotMat.at<Vec3b>(SRMrow, SRMcol)= SRMbgr;//assign color of pixel in SRM					
						} //done SRMcol
					} //done SRMrow
					//rotate the sun locus
					double xx = hcx + (scx - hcx)*cosgamma + (scy - hcy)*singamma;
					double yy = hcy - (scx - hcx)*singamma + (scy - hcy)*cosgamma;
					SunCenter.x = xx; SunCenter.y = yy;
					//crop to size with sun at center
					Rect roi;//region of interest centered at sun
					roi.x = int(SunCenter.x) - LSMCropSize/2.0; if (roi.x <0 ) {roi.x = 0;} if (roi.x >cols){ roi.x =cols; }//x offset
					roi.y = int(SunCenter.y) - LSMCropSize/2.0; if (roi.y < 0) {roi.y = 0;} if (roi.y >rows){ roi.y =rows; }//y offset
					roi.width = LSMCropSize; if ((roi.width + roi.x) > cols){ roi.width = cols-roi.x;} //width and height
					roi.height = LSMCropSize; if ((roi.height + roi.y) > rows){ roi.height = rows-roi.y;}
					Mat TempSkyMat = SkyRotMat(roi);	//create the cropped version, sun centered	but not past edge of image	
					//10-24-2018 at low sun angles, horizontal extension is stretched - I have not properly accounted for that.
					//ad hoc solution: scale the image by sun zenith angle
					int LSMxoffset=0, LSMyoffset = 0;
					if (ZenithAngle > MaximumZenithAngle/2)
					{
						double horizonscale = sin(MaximumZenithAngle/2.0*Pi/180.0)/sin(ZenithAngle*Pi/180);
						if (horizonscale < 0.7 ){ horizonscale = 0.7;}
						resize(TempSkyMat,TempSkyMat,Size(0,0),horizonscale,horizonscale); 
						LSMxoffset = 0.5*(1-horizonscale) * LSMCropSize; 
						LSMyoffset = 0.5*(1-horizonscale) * LSMCropSize; 
					}
					//redefine roi to center TempSkyMat in square LocalSkyMat
					roi.x = LSMxoffset; roi.y = LSMyoffset;
					roi.width = TempSkyMat.cols; roi.height = TempSkyMat.rows;
					TempSkyMat.copyTo(LocalSkyMat(roi)); //place the temporary crop into center of square LocalSkyMat
					SunCenter.x = LSMCropSize/2; SunCenter.y = LSMCropSize/2; // update the sun position as center of the LSM
					if (SaveLocalSkyMat == 1) //this is an option set in the header, avoid with big numbers of files
					{
						Mat SunTempMat; 
						LocalSkyMat.copyTo(SunTempMat);//initialize with the original image
						circle( SunTempMat, SunCenter, 6, Scalar( 0, 255, 255), -1, 10, 0); //the sun
						stringstream sc;sc.str(""); sc << CurrentImage<<"-7-LocalSkyMat.jpg"; string filename = sc.str(); imwrite(filename, SunTempMat);
						runlog << "\t saving LocalSkyMat in " << filename << "\n";
					}//end SaveLocalSkyMat
				} //****************** done create the Local Sky Map (LSM) ******************************************************

				
				//**************Collect the Average Radial Intensity around the sun in LocalSkyMat *******************************
				double ARI[3][4][3][200]; //matrix to hold the Average radial intensity distributions
				// index explanation:
				// 		[3]		color, 0=B, 1=G, 2=R
				//		[4]		quadrant, 0=TR, 1=BR, 2=BL, 3=TL
				//		[3]		type of data 0=ARI, 1=stdev(ARI), 2= number of points in quadrant
				//		[200]	ibin for solar distance in pixels. radius in sd = ibin * MaximumZenithAngle/HorizonRadius
				int SolDistMin = ShadowStripWidth / 2; // Minimum distance for radial analysis in pixel units
				int SolDistMax = LSMCropSize / 2; //Max distance for radial analysis in pixel units
				double ARIvar[3][4]; //get the average of the stdevs first index color, second index quadrant
				double ColorValue[4];// B2/RG average
				//			[4]		quadrant ( 0=TR, 1=BR, 2=BL, 3=TL )
				{// collect ARI
					double colortot[3][4];// collect sum of BGR for each color and quadrant
					for (int i=0; i<4; i++){ for (int j=0; j<3; j++) colortot[j][i] = 0.0;} //setting colortotal to zero
					for (int icolor = 0; icolor <3; icolor++)
					{

						for (int iquad=0; iquad <4; iquad++)
						{
							ARIvar[icolor][iquad] = 0.0;
							for (int idat =0; idat < 3; idat++)
							{
								for (int ibin = 0; ibin < 200; ibin++)
								{
									ARI[icolor][iquad][idat][ibin] = 0.0;
								}
							}
						}
					}
					double scx = SunCenter.x, scy = SunCenter.y; // pixel suncenter
					// scan all the pixels in the LocalSkyMat
					for (int irow = 0; irow < LocalSkyMat.rows; irow++ )
					{
						for (int icol = 0; icol < LocalSkyMat.cols; icol++)
						{
							double dx = double(icol) - scx, dy = double(irow) - scy;
							double radius = dx * dx + dy * dy; radius = sqrt(radius); // float pixel distance of local pixel from sun
							int SolDist = int(radius + 0.5); //solar distance of this pixel in integral pixel units
							if ((SolDist <= SolDistMax) && (SolDist >= SolDistMin))  //if within image and outsid of shadowstrip
							{//if within analysis area
								Vec3b bgrcolor = LocalSkyMat.at<Vec3b>(irow, icol);
								//assign the quadrant
								int iquad=0; //TR
								if ((dx > 0) && (dy > 0)) {iquad = 1;} //BR
								if ((dx < 0) && (dy > 0)) {iquad = 2;} //BL
								if ((dx < 0) && (dy < 0)) {iquad = 3;} //TL	
								for (int icolor=0; icolor<3; icolor++)
								{
									if ((icolor == 0 && bgrcolor[icolor] > 100) ||(icolor == 1 && bgrcolor[icolor] > 80) ||(icolor == 2 && bgrcolor[icolor] > 60)) //this limit should prevent remnant shadow strip from being counted?
									{
										//collect the colorvalue
										colortot[0][iquad] += bgrcolor[0];
										colortot[1][iquad] += bgrcolor[1];
										colortot[2][iquad] += bgrcolor[2];								
										for (int ibin = SolDist - SolDistBin; ibin < SolDist + SolDistBin; ibin++)
										{// all bins affecting this pixel
											if((SolDist >= (ibin - SolDistBin/2)) && (SolDist <= (ibin + SolDistBin/2)))
											{ //Assign to the correct radial distance
												ARI[icolor][iquad][0][ibin] +=  bgrcolor[icolor]; //intensity
												ARI[icolor][iquad][1][ibin] +=  bgrcolor[icolor]*bgrcolor[icolor]; //stdev
												ARI[icolor][iquad][2][ibin] +=  1.0; //counter
											} //done checking back and forth
										}// done all bins affecting this pixel
									} //done icolor conditions
								} //end icolor
							}//done within analysis area
						} //done icol						
					} //done irow
					//compute the color value B2/GR for each quadrant
					for (int iquad = 0; iquad < 4; iquad ++) 
					{ 
						ColorValue[iquad] = colortot[0][iquad] * colortot[0][iquad] / colortot[1][iquad] / colortot[2][iquad];
					} 	
					//compute the ARI and its properties
					for (int icolor = 0; icolor <3; icolor++)
					{
						for (int iquad=0; iquad <4; iquad++)
						{
							int ARIvarcount = 0;
							for (int ibin = 0; ibin <200; ibin++)
							{ //norm all bins
								double ave = 0.0; 
								if (ARI[icolor][iquad][2][ibin] != 0.0) {ave = ARI[icolor][iquad][0][ibin] / ARI[icolor][iquad][2][ibin];} //average intensity
								ARI[icolor][iquad][0][ibin] = ave; //assign average intensity
								double std = 0.0; 
								if (ARI[icolor][iquad][2][ibin] != 0.0) {std = sqrt(ARI[icolor][iquad][1][ibin] / ARI[icolor][iquad][2][ibin] - ave*ave);} //average intensity
								ARI[icolor][iquad][1][ibin] = std; //assign stdev
								if (std != 0.0)
								{
									ARIvar[icolor][iquad] += std;
									ARIvarcount += 1;
								}								
							} //done bins
							if (ARIvarcount != 0)
							{
								ARIvar[icolor][iquad] = ARIvar[icolor][iquad] / ARIvarcount;								
							}else
							{
								ARIvar[icolor][iquad] = 100.0;
							}
						} //done iquad
					} //done icolor
				}//************** Done Collect the Average Radial Intensity around the sun in LocalSkyMat *******************************
				
				
				//in order to eliminate problems of overexposure and wide shadow strips, etc, I define a FitStartRadius, separately for each color and each quadrant
				//the default value ofr FitStartRadius is rfitlimmin as read earlier from the parameter file
				double FitStartRadius [3][4];
				// indices: [3]	color	0-B, 1-G, 2=R-ARI
				//			[4]	quadrant	0-TR, 1-BR, 2-BL, 3-TL
				// now check for each color and quadrant, and assess the quality of the ARI data
					for (int icolor = 0; icolor < 3; icolor++)
					{
						for (int iquad = 0; iquad < 4; iquad++)
						{
							FitStartRadius[icolor][iquad] = rfitlimmin;
							int firstbin = int(rfitlimmin * HorizonRadius / MaximumZenithAngle + 0.5); // first bin
							int lastbin = int(rfitlimmax * HorizonRadius / MaximumZenithAngle + 0.5); // first bin
							if (ARI[icolor][iquad][0][firstbin] > 253.0 ) //likely overexposed
							{
								double intensityxx = 255.0;
								int ibin = firstbin;
								while ((intensityxx > 253.0) && (ibin < lastbin))
								{
									ibin+= 1;
									intensityxx = ARI[icolor][iquad][0][ibin];										
								}
								FitStartRadius[icolor][iquad] = double(ibin) * MaximumZenithAngle / HorizonRadius;
							}//done overexposed	
							double cutoffcolor = 100.0;
							if (icolor == 1) cutoffcolor = 80.0;
							if (icolor == 2) cutoffcolor = 50.0;
							if (ARI[icolor][iquad][0][firstbin] < cutoffcolor ) //likely parts of the shadow strip contained in ARI (occurs if not image not well-centered)
							{
								int pastmax = 0; //find first max and fit from there
								int ibin = firstbin;
								double intensityxx = 0.0;
								while ((( pastmax == 0 ) || ( intensityxx < 80.0 ) )|| (ibin < lastbin))
								{
									ibin+= 1;
									double intensityy = ARI[icolor][iquad][0][ibin];
									pastmax = 0;
                                    if ( intensityy > intensityxx )	{pastmax = 1;}
									intensityxx = intensityy;
								}
								FitStartRadius[icolor][iquad] = double(ibin) * MaximumZenithAngle / HorizonRadius;								
							}						
						} //done iquad
					}// done icolor
					
				//***** line fit of ARI for each quadrant and color *********************************************************************
				double LineFitPrms[3][4][3];
				// the index fields for LineFitPrms cover the following:
				//		[3]	colors (0-B, 1-G, 2-R)
				//		[4]	quadrant ( 0=TR, 1=BR, 2=BL, 3=TL )
				//		[3] line fit parameters 0 = slope, 1 = intercept, 2 = R2
				for (int icolor = 0; icolor < 3; icolor ++)
				{
					for (int iquad = 0; iquad < 4; iquad ++)
					{
						double slope, intercept, R2;
						// following reference mathworld.wolfram.com/LeastSquaresFitting.html
						double ssxx = 0.0, ssyy = 0.0, ssxy = 0.0; //variances and covariances of x and y values
						double xave = 0.0, yave = 0.0; //average values of x and y
						double xsd = FitStartRadius[icolor][iquad]; //starting radius in sky degrees
						int ibin = int(xsd * HorizonRadius / MaximumZenithAngle +0.5); //initial radial bin to evaluate (pixel value)
						double rfitnum = 0; //count the actual points used in fit
						while (xsd < rfitlimmax)
						{
							double x = xsd;
							double y = ARI[icolor][iquad][0][ibin];
							ssxx = ssxx + x*x;
							ssyy = ssyy + y*y;
							ssxy = ssxy + x*y;
							xave = xave + x;
							yave = yave + y;
							ibin = ibin + 1;
							xsd = double(ibin) * MaximumZenithAngle / HorizonRadius;
							rfitnum = rfitnum +1;
						}//done irad
						if (((iquad == 1 || iquad == 2 ) && ZenithAngle > 70.0 ) || rfitnum == 0 ) 
						{ // excludes low-sun situations
							slope = 1000.0; intercept = 1000.0; R2 =0.0;
						} else
						{
							xave = xave / rfitnum;
							yave = yave / rfitnum;
							ssxx = ssxx - rfitnum*xave*xave;
							ssyy = ssyy - rfitnum*yave*yave;
							ssxy = ssxy - rfitnum*xave*yave;
							slope = 0.0;
							if (ssxx != 0.0 ){slope = ssxy / ssxx;}							
							intercept = yave - slope*xave;
							R2 = ssxy*ssxy / ssxx / ssyy;	
						}
						LineFitPrms[icolor][iquad][0]=slope;
						LineFitPrms[icolor][iquad][1]=intercept;
						LineFitPrms[icolor][iquad][2]=R2;					
					}//done iquad
				}// done icolor
				if (SaveLineFitPRMs == 1)
				{
					linefitfile << CurrentImage <<", ";
					for (int iquad=0; iquad < 4; iquad++)
					{
						for (int icolor = 0; icolor < 3; icolor++)
						{
							for(int ifit =0; ifit < 3; ifit++)
							{
								linefitfile << LineFitPrms[icolor][iquad][ifit]<<",";
							}
						}
					}
					linefitfile << endl;
				}//done save line fit parameters
				//***** Done line fit of ARI for each quadrant and color *********************************************************************
				
				//*****************************************Analyze the running average of ARI, compute eta and etaderiv to characterize halo markers
				double ARI_RA[3][4][200]; //matrix to hold the running average of ARI
				int RALowBin = int(RALowLim * HorizonRadius / MaximumZenithAngle + 0.5);
				int RAHighBin = int(RAHighLim * HorizonRadius / MaximumZenithAngle + 0.5);
				double BinSize = MaximumZenithAngle / HorizonRadius; //angular size of one radial bin in sd
				double Eta[3][4][200]; //holds the data for ARI - ARI_RA
				double EtaDeriv[3][4][200]; //holds the values for the derivative dn/dr
				// index explanation:
				// 		[3]		color, 0=B, 1=G, 2=R
				//		[4]		quadrant, 0=TR, 1=BR, 2=BL, 3=TL
				//		[200]	ibin for solar distance in pixels. radius in sd = ibin * MaximumZenithAngle/HorizonRadius
				//compute the running average
				double RAwidthsd = 6.0; //defining the width over which the running average should be taken, 6 sd encompasses the whole halo width
				int RAwidthBins = int(RAwidthsd / BinSize +0.5); RAwidthBins = RAwidthBins / 2;  //ensures the range back and forth, resulting in even total.
				for (int icolor = 0; icolor < 3; icolor++)
				{
					for (int iquad =0; iquad < 4; iquad++ )
					{
						for (int ibin = SolDistMin + RAwidthBins; ibin <= SolDistMax - RAwidthBins; ibin++) //for all possible points
						{
							ARI_RA[icolor][iquad][ibin]=0.0;
							for (int i = -RAwidthBins; i <= RAwidthBins; i++)
							{
								ARI_RA[icolor][iquad][ibin] += ARI[icolor][iquad][0][ibin+i];
							}
							ARI_RA[icolor][iquad][ibin]=ARI_RA[icolor][iquad][ibin]/double(2* RAwidthBins +1);							
						}
						//get Eta = ARI - ARI_RA and EtaDeriv
						for (int ibin = RALowBin-1; ibin <= RAHighBin+1; ibin++)
						{
							Eta[icolor][iquad][ibin] = ARI[icolor][iquad][0][ibin] - ARI_RA[icolor][iquad][ibin];
						}
						for (int ibin = RALowBin; ibin <= RAHighBin; ibin++)
						{
							EtaDeriv[icolor][iquad][ibin] = ( Eta[icolor][iquad][ibin +1] - Eta[icolor][iquad][ibin - 1] )/BinSize/2.0 ;
						}
					}
				}
				if (SaveAveRadInt ==1 ) //the option is controlled in the header of the program, saves ARI, ARI_RA, and line in csv file
				{
					stringstream sc;sc.str(""); sc << CurrentImage<<"-AveRadInt.csv"; string filename = sc.str(); 
					ofstream ARIdata; ARIdata.open(filename);
					runlog << "\t saving Average Radial Intensities in " << filename << "\n";
					//write a header for the csv file
					ARIdata << " ** haloloop output for " << CurrentImage << " \n";
					ARIdata << " ** Average Radial Intensity (ARI) for each color channel versus radial distance from sun ***\n";
					ARIdata << " Minimum sampling radius:, " << SolDistMin << ",pixels,," << SolDistMin * MaximumZenithAngle/HorizonRadius << ", skydegrees \n";
					ARIdata << " Maximum sampling radius:, " << SolDistMax << ",pixels,," << SolDistMax * MaximumZenithAngle/HorizonRadius << ", skydegrees \n";
					ARIdata << " Running Average taken over , " << 2* RAwidthBins +1 << ",bins corresponding to approximately ," << RAwidthsd << ", skydegrees \n";
					ARIdata << " Sampling width:, " << SolDistBin << "\n four quadrants separately \n";
					ARIdata << " SolDist-pix, SolDist-sd,";
					ARIdata << " TR-B-ARI, TR-B-line, TR-B-RA, TR-B-eta, TR-B-dndr, TR-G-ARI, TR-G-line, TR-G-RA, TR-G-eta, TR-G-dndr,TR-R-ARI, TR-R-line, TR-R-RA, TR-R-eta, TR-R-dndr, ";  //done ARI
					ARIdata << " BR-B-ARI, BR-B-line, BR-B-RA, BR-B-eta, BR-B-dndr, BR-G-ARI, BR-G-line, BR-G-RA, BR-G-eta, BR-G-dndr,BR-R-ARI, BR-R-line, BR-R-RA, BR-R-eta, BR-R-dndr, ";  //done ARI
					ARIdata << " BL-B-ARI, BL-B-line, BL-B-RA, BL-B-eta, BL-B-dndr, BL-G-ARI, BL-G-line, BL-G-RA, BL-G-eta, BL-G-dndr,BL-R-ARI, BL-R-line, BL-R-RA, BL-R-eta, BL-R-dndr, ";  //done ARI
					ARIdata << " TL-B-ARI, TL-B-line, TL-B-RA, TL-B-eta, TL-B-dndr, TL-G-ARI, TL-G-line, TL-G-RA, TL-G-eta, TL-G-dndr,TL-R-ARI, TL-R-line, TL-R-RA, TL-R-eta, TL-R-dndr, ";  //done ARI
					ARIdata << " TR-B-stdev,  TR-G-stdev, TR-R-stdev, BR-B-stdev, BR-G-stdev, BR-R-stdev,  BL-B-stdev, BL-G-stdev, BL-R-stdev, TL-B-stdev, TL-G-stdev, TL-R-stdev, ";  //done ARI
					ARIdata << " TR-B-count,  TR-G-count, TR-R-count, BR-B-count, BR-G-count, BR-R-count,  BL-B-count, BL-G-count, BL-R-count, BL-T-count, TL-B-count, TL-G-count, TL-R-count \n";  //done ARI
					for (int ibin = SolDistMin; ibin <= SolDistMax; ibin ++ )
					{
						double radius = double(ibin) * MaximumZenithAngle / HorizonRadius;
						ARIdata << ibin << "," << radius << ",";             //corresponds to SolDist
						for (int iquad = 0; iquad < 4; iquad++)
						{
							for (int icolor=0; icolor<3; icolor++)
							{
								ARIdata << ARI[icolor][iquad][0][ibin] << ","; //ARI itself
								double linedat = LineFitPrms[icolor][iquad][0] * radius + LineFitPrms[icolor][iquad][1];
								ARIdata << linedat << ","; //the value of the local line fit
								ARIdata << ARI_RA[icolor][iquad][ibin] << ","; //the local running average RA
								ARIdata << Eta[icolor][iquad][ibin] << ","; // eta = ARI - RA
								ARIdata << EtaDeriv[icolor][iquad][ibin] << ","; //the derivative of eta dn/dr
							}//done icolor
						}//done iquad
						for (int idat=1; idat<3; idat++) //also save stdevs and counters as function of radius
						{
							for (int iquad = 0; iquad < 4; iquad++)
							{
								for (int icolor=0; icolor<3; icolor++)
								{
									ARIdata << ARI[icolor][iquad][idat][ibin] << ",";
								}//done icolor
							}//done iquad
						}//done idat
						ARIdata <<"\n";
					}//done ibin						
					ARIdata.close();
				} //done SaveAveRadInt	
						
				//********  Analysis id eta and eta derivative *********************************************
				double EtaDerivMarkers[3][4][9];
				// These markers will be used to assign the halo score after  this analysis segment
				//	index	def		meaning
				//	1		[3]		icolor (0-B, 1-G, 2-R)
				//	2		[4]		iquad (0-TR, 1-BR, 2-BL, 3-TL)
				//	3		[9]		Characteristics of the main maximum found in the ARI (largest change from etaderiv+ to etaderiv-)
				//					0	LocMaxDer
				//					1	LocPMDer
				//					2	LocMinDer
				//					3	ValMaxDer				
				//					4	ValMinDer			
				//					5	number of maxima detected in analysis interval				
				//					6	stdev across bgr for LocMaxDer (only icolor=0 slot used)
				//					7	stdev across bgr for LocPMDer (only icolor=0 slot used)
				//					8	stdev across bgr for LocMinDer (only icolor=0 slot used)				
				for (int iquad = 0; iquad < 4; iquad++)
				{
					for (int icolor = 0; icolor <3; icolor++)
					{
						double EtaDerivMaxProps[10][5];
						// holds the properties of maxima for the derivative of eta
						// 	index	def 	meaning
						//	1		[10]	index of the maxima found in the analysis interval, ideally one or two in halo pic
						//	2		[5]		a list of properties of this particular derivative maximum
						//					0	sd location of the maximum
						//					1	sd location of the subsequent pm crossing
						//					2	sd location of this subsequent minimum
						//					3	value of the maximum (in brightness change/sd)
						//					4	value of the subsequent minimum (brightness change/sd)
						int maxcounter = 0;
						double valmax = 0.0, maxloc = 0.0, pmloc = 0.0, minloc = 0.0, valmin = 0.0;
						//find the first positive value
						int firstposbin = RALowBin;
						while (EtaDeriv[icolor][iquad][firstposbin] < 0 ){ firstposbin += 1;}
						int ibin = firstposbin;
						int posrange = 1; //will swith to -1 if pm cross
						while (ibin < RAHighBin)
						{
							valmax = EtaDeriv[icolor][iquad][ibin]; 
							maxloc = double(ibin);
							int maxfound =0;
							while (posrange == 1) // move forward through the entire positive range
							{
								double currentdndr = EtaDeriv[icolor][iquad][ibin];							
								if (currentdndr > valmax) 
								{	//collect the max value and location in the positive range
									valmax = currentdndr;
									maxloc = double(ibin);
								}
								if (currentdndr < 0.0)
								{	//positive range ended, so reset and save what you know
									posrange = -1;
									maxfound =1;
									pmloc =double(ibin)-0.5;
									maxcounter += 1;												
									EtaDerivMaxProps[maxcounter][3] = valmax;
									EtaDerivMaxProps[maxcounter][0] = maxloc * MaximumZenithAngle / HorizonRadius;
									EtaDerivMaxProps[maxcounter][1] = pmloc * MaximumZenithAngle / HorizonRadius;
								}
								ibin += 1;
								if (ibin == RAHighBin) posrange = 0;
							}//end of positive range
							ibin = ibin - 1;
							if (maxfound == 1) {valmin = 0.0;	}
							int minfound = 0;														
							while ((minfound == 0)&& (ibin <= RAHighBin))//find the first minimum
							{
								double currentdndr = EtaDeriv[icolor][iquad][ibin];															
								if ( currentdndr < valmin ) 
								{
									valmin = currentdndr;
									minloc = double(ibin);									
								} 
								if ( currentdndr > (0.8 * valmin )) //the offset may prevent local noise from registering an early minimum
								{
									minfound = 1;
								}									
								ibin += 1;
							}
							//save the location and value for this minimum
							if ( minfound == 1 )
							{
								EtaDerivMaxProps[maxcounter][4] =  valmin;
								EtaDerivMaxProps[maxcounter][2] = minloc * MaximumZenithAngle / HorizonRadius;
							}
							//sometimes, no subsequent minimum is found. In that case LocMinDer and ValMinDer need special handling
							if ((minfound == 0)&&(maxfound==1))
							{
								EtaDerivMaxProps[maxcounter][4] =  EtaDeriv[icolor][iquad][RAHighBin];
								EtaDerivMaxProps[maxcounter][2] = RAHighBin * MaximumZenithAngle / HorizonRadius;
							}
							while ( posrange == -1 ) 
							{
								if (EtaDeriv[icolor][iquad][ibin] > 0.0) posrange = 1;
								if (ibin > RAHighBin)  posrange = 0;
								ibin +=1 ;// now move forward through the entire negative range
							}
						}//while in analysis interval (ibin < RAHighBin)
						
						// //give me some output
						// for (int imax = 1; imax <= maxcounter; imax++)
						// {
							// cout << "image, icolor, iquad, imax,  maxloc, pmloc, minloc, valmax, valmin:  "<< CurrentImage << "\t" << icolor << "\t" << iquad << "\t" << imax;
							// for (int idat = 0; idat < 5; idat++)
							// {
								// cout << " \t" << EtaDerivMaxProps[imax][idat];
							// }
							// cout << "\n";
						// }
						// cout << "\n";
						
						//decide on the main maximum of etaderiv found and populate the field EtaDerivMarkers [0-2]
						double LocMaxDer = 0.0, LocPMDer = 0.0, LocMinDer = 0.0, ValMaxDer = 0.0, ValMinDer = 0.0;
						int maxrangeindex = 0;
						double maxrange=0.0;
						for (int imax = 1; imax <= maxcounter; imax++) //find the indexed maximum with the largest index
						{
							double range = EtaDerivMaxProps[imax][3] - EtaDerivMaxProps[imax][4];
							if (range > maxrange) 
							{
								maxrange = range;
								maxrangeindex = imax;
								LocMaxDer = EtaDerivMaxProps[maxrangeindex][0];
								LocPMDer = EtaDerivMaxProps[maxrangeindex][1];
								LocMinDer = EtaDerivMaxProps[maxrangeindex][2];
								ValMaxDer = EtaDerivMaxProps[maxrangeindex][3];
								ValMinDer = EtaDerivMaxProps[maxrangeindex][4];
							}
						}
						//cout<< " etaderivmarkers \t "  << "\t" << icolor << "\t" << iquad << "\t" << LocMaxDer << "\t" << LocPMDer << " \t" << LocMinDer << "\n";
						EtaDerivMarkers[icolor][iquad][0] = LocMaxDer; 	// sd position of the maximum derivative/ sd position of the upslope toward halo
						EtaDerivMarkers[icolor][iquad][1] = LocPMDer; 	// sd position of the zero derivative / sd position of max halo brightness
						EtaDerivMarkers[icolor][iquad][2] = LocMinDer; 	// sd position of the maximum negative derivative / sd position of downslope after halo						
						EtaDerivMarkers[icolor][iquad][3] = ValMaxDer; //number of maxima found, in halo image should be small						
						EtaDerivMarkers[icolor][iquad][4] = ValMinDer; //number of maxima found, in halo image should be small						
						EtaDerivMarkers[icolor][iquad][5] = double(maxcounter); //number of maxima found, in halo image should be small						
					}//done icolor
					//correlate the markers; in a halo image all markers are consistent within bin resolution with each other
					for (int imarker = 0; imarker < 3; imarker++ )
					{
						double sum = 0.0, sqsum = 0.0, stdev = 0.0;
						double counter=0;
						for (int icolor = 0; icolor < 3; icolor++)
						{
							double value = EtaDerivMarkers[icolor][iquad][imarker];
							if (value != 0.0 )
							{
								sum += value; 
								sqsum += value * value;		
								counter += 1.0;
							}
						}
						if (counter > 0.0 )
						{
							stdev = (sqsum / counter) - (sum * sum / counter / counter); 
							if (stdev > 0.0) {stdev = sqrt(stdev);} 	
						}
						if (stdev < 1.e-4){ stdev = 0.0;}
						EtaDerivMarkers[0][iquad][imarker + 6] = stdev;
						EtaDerivMarkers[1][iquad][imarker + 6] = stdev;
						EtaDerivMarkers[2][iquad][imarker + 6] = stdev;
					}
				}// done iquad
				// Done  ********  Analysis id eta and eta derivative *********************************************
				
				//****  Save the property file (lines can be used to improve the scoring by copy/past to Haloloop-Skytype-HaloScore-Parametercollection.xlsx )
				if (SavePropertyLog == 1)
				{
					for (int iquad = 0; iquad < 4; iquad++)
					{
						PropertyLog << CurrentImage << "," << iquad << "," ;
						for (int icolor =0; icolor <3; icolor++ )
						{
							PropertyLog << LineFitPrms[icolor][iquad][0] <<  ","; 	//slope
							PropertyLog << LineFitPrms[icolor][iquad][1] << ",";		//intercept
							PropertyLog << ARIvar[icolor][iquad] << ",";				//ARIvar
							PropertyLog << ColorValue[iquad] << ","; 				//CV
							for (int imarker = 0; imarker <9; imarker++)
							{
								PropertyLog << EtaDerivMarkers[icolor][iquad][imarker]<< ",";		//LocMaxDer, LocPMDer, LocMinDer, ValMaxDer, ValMinDer, nmax, stdLMaD, stdLPMD, stdLMiD							
							}
						}//done icolor
						for (int icolor =0; icolor <3; icolor++ )
						{
							PropertyLog << EtaDerivMarkers[icolor][iquad][1] - EtaDerivMarkers[icolor][iquad][0] << ","; //mxtozero
							PropertyLog << EtaDerivMarkers[icolor][iquad][2] - EtaDerivMarkers[icolor][iquad][0] << ","; //maxtomin							
						}
						for (int icolor =0; icolor <3; icolor++ )
						{
							PropertyLog << EtaDerivMarkers[icolor][iquad][3] - EtaDerivMarkers[icolor][iquad][4] << ","; //DerRange
						}
						PropertyLog << "\n";
					}//done iquad					
				}// *** done save property log
				
				//**** form the property vectors X for both, skytype and haloscore
				double SkyTypePrmsX[4][10]; // the vector of the current values for properties relevant in  sky type scoring
				//	index		meaning
				//	1	[4]		iquad (0-TR, 1-BR, 2-BL, 3-TL)
				//	2	[10]	property index:
				//					0-	B	 slope
				//					1-	B	 intercept
				//					2-	B	 ARIvar
				//					3-	B	 CV
				//					4-	G	 slope
				//					5-	G	 intercept
				//					6-	G	 ARIvar
				//					7-	R	 slope
				//					8-	R	 intercept
				//					9-	R	 ARIvar
				double HaloScorePrmsX[4][31]; // the vector of current property values relevant for the halo score
				//	index		meaning
				//	1	[4]		iquad (0-TR, 1-BR, 2-BL, 3-TL)
				//	2	[31]	property index:
				//				0	B slope				1	B intercept			2	B ARIvar			3	CV				4	B LocMaxDer
				//				5	B LocPMDer			6	B LocMinDer			7	B ValMaxDer			8	B ValMinDer		9	B nmax
				//				10	B stdLMaD			11	B stdLPMD			12	B stdLMiD			13	G slope			14	G intercept
				//				15	G ARIvar			16	G LocMaxDer			17	G LocPMDer			18	G LocMinDer		19	G ValMaxDer
				//				20	G ValMinDer			21	G nmax				22	R slope				23	R intercept		24	R ARIvar
				//				25	R LocMaxDer			26	R LocPMDer			27	R LocMinDer			28	R ValMaxDer		29	R ValMinDer
				//				30	R nmax							
				for (int iquad = 0; iquad < 4; iquad++)
				{
					SkyTypePrmsX[iquad][0] = LineFitPrms[0][iquad][0]; 	// B slope
					SkyTypePrmsX[iquad][1] = LineFitPrms[0][iquad][1]; 	// B intercept
					SkyTypePrmsX[iquad][2] = ARIvar[0][iquad];			// B ARIvar
					SkyTypePrmsX[iquad][3] = ColorValue[iquad];			// CV
					SkyTypePrmsX[iquad][4] = LineFitPrms[1][iquad][0]; 	// G slope
					SkyTypePrmsX[iquad][5] = LineFitPrms[1][iquad][1]; 	// G intercept
					SkyTypePrmsX[iquad][6] = ARIvar[1][iquad];			// G ARIvar
					SkyTypePrmsX[iquad][7] = LineFitPrms[2][iquad][0]; 	// R slope
					SkyTypePrmsX[iquad][8] = LineFitPrms[2][iquad][1]; 	// R intercept
					SkyTypePrmsX[iquad][9] = ARIvar[2][iquad];			// R ARIvar					
					//for (int iprop = 0; iprop <10; iprop ++) { cout << SkyTypePrmsX[iquad][iprop] << "\t";}
					//cout << "\n";
					HaloScorePrmsX[iquad][0] = LineFitPrms[0][iquad][0]; 		// B slope
					HaloScorePrmsX[iquad][1] = LineFitPrms[0][iquad][1]; 		// B intercept
					HaloScorePrmsX[iquad][2] = ARIvar[0][iquad];				// B ARIvar
					HaloScorePrmsX[iquad][3] = ColorValue[iquad];				// CV
					HaloScorePrmsX[iquad][4] = EtaDerivMarkers[0][iquad][0];	// B LocMaxDer
					HaloScorePrmsX[iquad][5] = EtaDerivMarkers[0][iquad][1];	// B LocPMDer
					HaloScorePrmsX[iquad][6] = EtaDerivMarkers[0][iquad][2];	// B LocMinDer
					HaloScorePrmsX[iquad][7] = EtaDerivMarkers[0][iquad][3];	// B ValMaxDer
					HaloScorePrmsX[iquad][8] = EtaDerivMarkers[0][iquad][4];	// B ValMinDer
					HaloScorePrmsX[iquad][9] = EtaDerivMarkers[0][iquad][5];	// B nmax
					HaloScorePrmsX[iquad][10] = EtaDerivMarkers[0][iquad][6];	// stdLMaD
					HaloScorePrmsX[iquad][11] = EtaDerivMarkers[0][iquad][7];	// stdLPMD
					HaloScorePrmsX[iquad][12] = EtaDerivMarkers[0][iquad][8];	// stdLMiD
					HaloScorePrmsX[iquad][13] = LineFitPrms[1][iquad][0]; 		// G slope
					HaloScorePrmsX[iquad][14] = LineFitPrms[1][iquad][1]; 		// G intercept
					HaloScorePrmsX[iquad][15] = ARIvar[1][iquad];				// G ARIvar
					HaloScorePrmsX[iquad][16] = EtaDerivMarkers[1][iquad][0];	// G LocMaxDer
					HaloScorePrmsX[iquad][17] = EtaDerivMarkers[1][iquad][1];	// G LocPMDer
					HaloScorePrmsX[iquad][18] = EtaDerivMarkers[1][iquad][2];	// G LocMinDer
					HaloScorePrmsX[iquad][19] = EtaDerivMarkers[1][iquad][3];	// G ValMaxDer
					HaloScorePrmsX[iquad][20] = EtaDerivMarkers[1][iquad][4];	// G ValMinDer
					HaloScorePrmsX[iquad][21] = EtaDerivMarkers[1][iquad][5];	// G nmax
					HaloScorePrmsX[iquad][22] = LineFitPrms[2][iquad][0]; 		// R slope
					HaloScorePrmsX[iquad][23] = LineFitPrms[2][iquad][1]; 		// R intercept
					HaloScorePrmsX[iquad][24] = ARIvar[2][iquad];				// R ARIvar
					HaloScorePrmsX[iquad][25] = EtaDerivMarkers[2][iquad][0];	// R LocMaxDer
					HaloScorePrmsX[iquad][26] = EtaDerivMarkers[2][iquad][1];	// R LocPMDer
					HaloScorePrmsX[iquad][27] = EtaDerivMarkers[2][iquad][2];	// R LocMinDer
					HaloScorePrmsX[iquad][28] = EtaDerivMarkers[2][iquad][3];	// R ValMaxDer
					HaloScorePrmsX[iquad][29] = EtaDerivMarkers[2][iquad][4];	// R ValMinDer
					HaloScorePrmsX[iquad][30] = EtaDerivMarkers[2][iquad][5];	// R nmax
				}//done iquad

				//****************** Compute the SkyType Scores! ******************************************
				double SkyTypeScoreRaw[4][4][2]; // stores the raw scores after matrix manipulation and normal distribution
				// index explanation
				//	1	[4]		iquad (0-TR, 1-BR, 2-BL, 3-TL)
				//	2	[4]		sky type index: 0-CS, 1-PCL, 2-CLD, 3-CLR, 
				//	3	[2]		0 - exponential, 1 - Mahalanobis distance squared (just for checks and balances)
				double SSTfactor = 1.0e3; // the prefactor for the Gaussian, it's arbitrary but tested in spreasdheet to give some doable numbers
				double SkyTypeScoreRel[4][5]; // stores the relative skyscores, as the quadrant score sum is set to 100%
				// index explanation
				//	1	[4]		iquad (0-TR, 1-BR, 2-BL, 3-TL)
				//	2	[5]		sky type index: 0-CS, 1-PCL, 2-CLD, 3-CLR, 4-none/not determined
				string DominantSky[5]; DominantSky[0] = "CS "; DominantSky[1] = "PCL"; DominantSky[2] = "CLD"; DominantSky[3] = "CLR"; DominantSky[4] = "---";
				int dominantsky[4];//save the domminant itype for each quadrant				
				for (int iquad =0; iquad < 4; iquad++ )
				{
					for (int itype = 0; itype < 4; itype++)
					{
						double XM[10]; //(X-M)
						for (int iprop = 0; iprop < 10; iprop++) //construct (X-M)
						{
							XM[iprop] = SkyTypePrmsX[iquad][iprop]-SkyTypePrmsM[itype][iprop];
						}
						double SXM[10];//S^(-1)(X-M) part for skytype Mahalanobis distance 
						for (int iprop =0; iprop < 10; iprop++ ) //construct S^(-1)(X-M)
						{
							SXM[iprop] =0;
							for (int jprop = 0; jprop <10; jprop++)
							{
								SXM[iprop]+= SkyTypePrmsICM[itype][iprop][jprop]*XM[jprop];//S^(-1)(X-M)
							}//done summing over iprop							
						}//done jprop
						double D2 = 0.0; //Mahalanobis distance D2 = (X-M)^T S^(-1)(X-M)
						for (int iprop =0; iprop < 10; iprop++ )
						{					
							D2+= XM[iprop]* SXM[iprop];
						}
						SkyTypeScoreRaw[iquad][itype][0] = SSTfactor * exp(-0.5*D2);
						SkyTypeScoreRaw[iquad][itype][1] = D2;
					}//done itype
					//get the relative sky scores for this quadrant
					double scoresum = 0.0;
					double maxscore = 0.0;
					dominantsky[iquad]=4;//default setting to N/A
					SkyTypeScoreRel[iquad][4] = 100.0;//default setting for undetermined skytype
					for (int itype = 0; itype < 4; itype++ )
					{
						double xx = SkyTypeScoreRaw[iquad][itype][0];
						scoresum += xx;
						if (xx > maxscore) 
						{
							maxscore = xx;
							dominantsky[iquad] = itype;
						}
					}
					if ((scoresum > 1.e-8) && (dominantsky[iquad] <4))
					{
						for (int itype = 0; itype < 4; itype++ )
						{
							SkyTypeScoreRel[iquad][itype] = SkyTypeScoreRaw[iquad][itype][0] / scoresum * 100.0;
						}							
						SkyTypeScoreRel[iquad][4] = 0.0;//reset undetermined skytype if an actual skytype has been found
					} else
					{
						for (int itype = 0; itype < 4; itype++ )
						{
							SkyTypeScoreRel[iquad][itype] = 0.0;
							dominantsky[iquad] = 4;
						}													
					}
					// give me some output
					// cout << CurrentImage << "\t" << iquad << "\t" << DominantSky[dominantsky[iquad]] << "\n ";
					// for (int itype = 0; itype < 5; itype++ )
					// {
						// cout << "\t \t \t " << SkyTypeScoreRaw[iquad][itype][0] << "\t" << SkyTypeScoreRel[iquad][itype] << " % for " << DominantSky[itype] <<" \t \n";
					// }
					// cout << "\n";
				}//done iquad	
				//done ****************** Compute the SkyType Scores! ******************************************
				
				
				//****************** Compute the Halo Scores! ******************************************
				double HaloScore[4][2]; // stores the halo score after matrix manipulation 
				// index explanation
				//	1	[4]		iquad (0-TR, 1-BR, 2-BL, 3-TL)
				// 	2	[2]		0=haloscore, 1=D Mahalanobis distance
				for (int iquad =0; iquad < 4; iquad++ )
				{
					double XM[31]; //(X-M)
					for (int iprop = 0; iprop < 31; iprop++) //construct (X-M)
					{
						XM[iprop] = HaloScorePrmsX[iquad][iprop]-HaloScorePrmsM[iprop];
					}
					double SXM[31];//S^(-1)(X-M) part for skytype Mahalanobis distance 
					for (int iprop =0; iprop < 31; iprop++ ) //construct S^(-1)(X-M)
					{
						SXM[iprop] =0.0;
						for (int jprop = 0; jprop <31; jprop++)
						{
							SXM[iprop]+= HaloScorePrmsICM[iprop][jprop]*XM[jprop];//S^(-1)(X-M)
						}//done summing over iprop							
					}//done jprop
					double D2 = 0.0; //Mahalanobis distance D2 = (X-M)^T S^(-1)(X-M)
					for (int iprop =0; iprop < 31; iprop++ )
					{					
						D2+= XM[iprop]* SXM[iprop];
					}
					if (D2 < 0) { D2 = 500.0; }
					if (D2 > 1000.0) { D2 = 500.0; }
					
					//score the D2, based on measurements taken in I:\haloloop\test\parameter-testing\RAtesting\haloloop-out-PropertyLog-sorting.xlsx for HaloScore
					HaloScore[iquad][0] = 0.0;
					if (D2 < 100.0){HaloScore[iquad][0] = 1.e6 * exp(-0.5*D2) ;} //compute the Gaussian
					HaloScore[iquad][1] = D2 ; //compute the Gaussian
					// give me some output
					//cout << CurrentImage << "\t" << iquad << "\t" << HaloScore[iquad][0] << "\t" << HaloScore[iquad][1]<< " \n";
				}//done iquad	
				//done ****************** Compute the Halo Scores! ******************************************

				//Now handle the output to SkyTypeLog
				//assign a general skytype:
				int iGenST = 4; // setting default for General Skytype to 4 = N/A
				double sumsky[5]; 
				for (int itype = 0; itype < 5; itype++){sumsky[itype]=0.0;} // sum for all sky types
				for (int iquad =0; iquad < 4; iquad++)
				{
					for (int itype = 0; itype < 5; itype++) // including N/A
					{
						sumsky[itype] += SkyTypeScoreRel[iquad][itype];
					}
				}
				for (int itype = 0; itype < 5; itype++)
				{
					sumsky[itype] = sumsky[itype]/4.0; //averageing over all quadrants
				}
				//gather the max
				double skymax = 0.0;
				for (int itype = 0; itype < 5; itype++)
				{
					if (sumsky[itype] > skymax)
					{
						skymax = sumsky[itype];
						iGenST = itype;
					}
					//cout << "\t total \t " << sumsky[itype] << " % for " << DominantSky[itype] <<" \t \n";
				}

				//assign a total haloscore as the sum over all quadrants (it's not a percentage!)
				double halosum = 0.0;
				for (int iquad = 0; iquad <4; iquad++)
				{
					halosum += HaloScore[iquad][0];
				}
				string Halo = "---";
				if (halosum > 1.0) {Halo = "HAL";}
				if ((halosum > 1.e-4 ) && (halosum < 1.0)) {Halo = "hal";}
				
				
				//compute the current time to handle gaps in the recorded, in h
				double CurrentTime = double(day)*24.0 + double(hour) +double(minute)/60.0 + double(second)/3600.0;
				if (imagecounter == 0 ) {PreviousTime = CurrentTime; }
				if ((CurrentTime - PreviousTime) > 1.0) //a new time series is starting: empty line, flush buffer
				{
					// flush the buffer and restart it for a new time series
					int itop= 100 + ContImage;
					if (ContImage > BuffRange) itop = 100 + BuffRange;
					for (int iwrite = itop; iwrite > 100; iwrite-- )
					{
						SkyTypeFile << BufferString[iwrite][0] << ",";
						for (int iquad = 0; iquad <5; iquad ++ ) { SkyTypeFile << BufferDos [iwrite][iquad]<< ",";}					
						SkyTypeFile << BufferString[iwrite][1] << "\n";
					}
					for (int ibuff = 100 - BuffRange; ibuff <= 100 + BuffRange; ibuff++)
					{
						for (int iquad = 0; iquad < 5; iquad++)
						{
							BufferDos[ibuff][iquad]=0.0;
						}
					}
					SkyTypeFile << "\n";
					HaloLogFile << "\n";
					ContImage = 0;					
				}
				PreviousTime = CurrentTime; //reset the time comparison
				ContImage = ContImage +1; // move continuous image counter forward. always 1 for the first i,age of a series
			
				//collect the findings into a buffer to enable a dos computation. That is basically a broadening of every haloscore across subsequent images to combine them into a halodos.
					// assess where buffer bottom and buffer top are. Higher buffer indexes for older images
					BuffTop = 100 + ContImage-1;//buffer is not filled to 3sigma range yet
					if (ContImage > BuffRange) BuffTop = 100 + BuffRange; // 
					BuffBot = 100 - BuffRange;
					
				
					//create the string info for the current image, current image is always held at buffer position 100
					sc << origin << "," << month << "/"  << day << "/" << year   << " "<< hour << ":" << minute << ":" << second << "," ;
					sc << CurrentImage << "," << AzimuthAngle << "," << ZenithAngle << ",," <<  DominantSky[iGenST] << "," << Halo << ",";
					BufferString[100][0] = sc.str(); sc.str(""); //skytypelog line before halodos				
					sc <<"," << sumsky[0] << ","<< sumsky[1] << ","<< sumsky[2] << ","<< sumsky[3] << ","<< sumsky[4] <<  "," << halosum << ",,";				
					for (int iquad = 0; iquad <4; iquad++)
					{
						sc << DominantSky[dominantsky[iquad]] << ",";
						for (int itype = 0; itype < 5; itype++ )
						{
							sc << SkyTypeScoreRel[iquad][itype] << ",";
						}
						sc << HaloScore[iquad][0]<< "," << HaloScore[iquad][1]<< ",,";
					}
					BufferString[100][1] = sc.str(); sc.str(""); //skytypelog line after halodos
			
					//create the Gaussian for the current haloscore and add it to the buffer
					for (int ibuff = BuffBot; ibuff <= BuffTop; ibuff++ )
					{
						double tdelay = 100.0 - double(ibuff);
						for (int iquad = 0; iquad <4; iquad++)
						{
							BufferDos[ibuff][iquad] += HaloScore[iquad][0]* exp(-0.5 * tdelay * tdelay / HalfWidth / HalfWidth); // halo dos for each quadrant
						}
						BufferDos[ibuff][4]+= halosum * exp(-0.5 * tdelay * tdelay / HalfWidth / HalfWidth); // all sky halo dos				
					}
			
					// write some of the buffer out to the respective log files.
					if (ContImage >= BuffRange )//image should have progressed to top of buffer before written out
					{
						SkyTypeFile << BufferString[BuffTop][0] << ",";
						for (int iquad = 0; iquad <5; iquad ++ ) { SkyTypeFile << BufferDos [BuffTop][iquad]<< ",";}					
						SkyTypeFile << BufferString[BuffTop][1] << "\n";
					}
				
					// push the buffer entries forward
					for (int ibuff = 100 + BuffRange; ibuff > 100 - BuffRange ; ibuff-- )
					{
						BufferString[ibuff][0] = BufferString[ibuff-1][0];
						BufferString[ibuff][1] = BufferString[ibuff-1][1];
						for (int j=0; j<5; j++) {BufferDos[ibuff][j] = BufferDos[ibuff-1][j];}
					}
					for (int j=0; j<5; j++) {BufferDos[100-BuffRange][j] =0.00;} //fill first line with zero
			
				//make a note in runlog
				runlog << " \t all-sky type: \t" << DominantSky[iGenST] << " \t All-sky halo: \t" <<Halo << "\n";
	
				
			}// end else if (!SourceImageMat.data) = image was not empty
		} //end for (imagecounter = 0; imagecounter < totalimagecount; imagecounter++) 
	//******************************End processing images *********************************************************
	
	//Flush the buffer
	for (int ibuff = 100 + BuffRange; ibuff > 100; ibuff--)
	{
		SkyTypeFile << BufferString[ibuff][0] << ",";
		for (int iquad = 0; iquad <5; iquad ++ ) { SkyTypeFile << BufferDos [ibuff][iquad]<< ",";}					
		SkyTypeFile << BufferString[ibuff][1] << "\n";
	}
	
	
	
	//***********************  mopping up *********************************************
	clock_t endTime = clock();
	clock_t clockticks = endTime - startTime;
	double time = clockticks / CLOCKS_PER_SEC;
	double timemin = time / 60.0;
	cout << "\n\nTotal number of images processed: " << imagecounter << endl;
	cout << "Time was " << timemin << " minutes." << endl;
	runlog << "\n\nTotal number of images processed: " << imagecounter << endl;
	runlog << "Time was " << timemin << " minutes." << endl;
	if (SaveLineFitPRMs == 1) linefitfile.close();
	if (SavePropertyLog == 1) PropertyLog.close();
	SkyTypeFile.close();
	runlog.close();

	
}