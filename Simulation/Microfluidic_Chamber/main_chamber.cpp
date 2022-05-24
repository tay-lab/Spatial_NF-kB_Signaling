/****************************************
* Main file executing the simulation    *
****************************************/
/* This program generates data from a   *
simulation of NFkB activation. The acti-*
vation is caused by a diffusing concen- *
tration of TNF, which is produced by    *
different types of sources.             *
****************************************/

// Include:
#include "cell_chamber.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "MersenneTwister.h"	
#include <assert.h>
#include <time.h>
#include "omp.h"
using namespace std;

// Main function:
int main(){

	/* Start timer. The final run time is displayed in the command window */
	double start = omp_get_wtime();
	
	Cell C[Xcells][Ycells];			// Ncells number of objects of type cell. Cell contains all differential equations and the RK solver.
	double Tnew[Xcells][Ycells];	// Matrix containing new cytokine values after diffusion.
	double d = 1200.;				// Diffusion constant.
	double taui = 45./60.;			// Half-life TNF in hours.
	double newS1 = 7.5; 			// Strength of signal 1.
	string nS = "7_5";
	int preChamber = 30;			// Size of prechamber.
	string ens = "10";				// Name of the given simulation.
	
	MTRand(time(NULL));	// Seed the random number generator. The random number generator is the Mersenne Twister.
	MTRand RN;			// Create random number as a class object of class MTRand.

	// Define time:
	double tin, tstart, tfin, dt, t, tstop;
	tin = -10; tfin = 10; 				// Initial and final time of the simulations.
	tstart = -0.0001; tstop = 0.0001; 	// 0.0001	// Time to turn source on and off.
	dt = 0.0001;
	int PrintDummy = 0.05/dt; 			// PrintDummy is the frequency of saving data to files.
	string nT = "dt";					// Duration of TNF stimulus
	double RanNum;
	
	// Assign initial values:
	for (int jCell = 0; jCell<Ycells; jCell++){
		for (int iCell = 0; iCell<Xcells; iCell++){
			C[iCell][jCell].AssignParameters();															// Assign parameters to each cell.
			C[iCell][jCell].Initialize(0.0436,    0.1124,    4.0468,    0.0095,    1.3770, 0., 0.); 	// Initialize all variables in each cell. 
			// initNF, initIm, initIa, initA20m, initA20, initIKK, initIKKi, determined in Matlab by running for 300 hr until equilibrated.
			RanNum = RN.randNorm(1.0,0.2);
			if (RanNum<0){RanNum=0.001;}		  // Negative values are set to 0.001.
			C[iCell][jCell].ReAssignKTNF(RanNum); // Each cell is assigned a random activation threshold.
		}
	}
	
/*******************************************************************/	
	
	// Generate txt files for saving data:
	ofstream myfile;
	ofstream myfile1;
	ofstream myfile2;
	ofstream myfile3;
	ofstream myfile4;
	ofstream myfile5;
	
	string fileName1, fileName2, GName1, GName2, GName3, fileName3, filepath;
	fileName1 = "NFkB_time_" + nT + "_S_" + nS + "_" + ens +".txt";
	fileName2 = "TNF_time_" + nT + "_S_" + nS + "_" + ens +".txt";
	GName1 = "Gene1_time_"+nT+"_S_"+nS+"_" + ens +".txt";
	GName2 = "Gene2_time_"+nT+"_S_"+nS+"_" + ens +".txt";
	GName3 = "Gene3_time_"+nT+"_S_"+nS+"_" + ens +".txt";
	fileName3 = "TotTNF_time_" +nT+ "_S_" +nS+ "_" + ens +".txt"; 
	filepath = "DATA\\";
	
	myfile.open((filepath+fileName1).c_str());
	myfile1.open((filepath+fileName2).c_str());
	myfile2.open((filepath+GName1).c_str());
	myfile3.open((filepath+GName2).c_str());
	myfile4.open((filepath+GName3).c_str());
	myfile5.open((filepath+fileName3).c_str());

/**********************************************************
*** Perform timesteps, this is where the program starts: **
***********************************************************/
	
	t = tin;
	double TNFhandle = 0;
	
for (int i=0;t<tfin;i++){
	
	// Reset gene expression so gene expression due to the equilibration process is not taken into account.
	if (abs(t-tstart)<1e-5){
		for (int jCell = 0; jCell<Ycells; jCell++){
			for (int iCell = 0; iCell<Xcells; iCell++){
				C[iCell][jCell].InitializeGenes(0.0001);
			}
		}
	}

	// Turn on source if t is between tstart and tstop:
	if (((t-tstart)>0) && ((t-tstop)<0)){
		for (int iCell = 0; iCell < Xcells; iCell++){
			for (int jCell = 0; jCell < preChamber; jCell++){
				C[iCell][jCell].ReassignT(newS1);	// The prechamber is filled with TNF.
			}										
		}
	}

	// Send TNF values of each cell to the txt.file.
	if (i%PrintDummy == 0){
	TNFhandle = 0;	
		
	myfile  << t << " ";
	myfile1 << t << " ";
	myfile2 << t << " ";
	myfile3 << t << " ";
	myfile4 << t << " ";
	myfile5 << t << " ";
		
	for (int jCell = 0; jCell<Ycells; jCell++){
		for (int iCell = 0; iCell<Xcells; iCell++){
			myfile  << C[iCell][jCell].NOut()  << " ";
			myfile1 << C[iCell][jCell].TOut()  << " ";
			myfile2 << C[iCell][jCell].G1Out() << " ";
			myfile3 << C[iCell][jCell].G2Out() << " ";
			myfile4 << C[iCell][jCell].G3Out() << " ";
			TNFhandle = TNFhandle + C[iCell][jCell].TOut();
		}
	};
	
	myfile  << endl;
	myfile1 << endl;
	myfile2 << endl;
	myfile3 << endl;
	myfile4 << endl;
	myfile5 << TNFhandle/((Xcells)*(Ycells)) << endl;
	}
	
	// Update all the cells:
	int iCell;
	# pragma omp parallel for private (iCell)	// Do updates in parallel.
	for (int jCell = preChamber; jCell<(Ycells); jCell++){
		for (int iCell = 0; iCell<(Xcells); iCell++){
				C[iCell][jCell].Timestep(dt);
		}
	}

	// Diffusion and decay of T is calculated:
	int jCell;
	# pragma omp parallel for private (jCell) 	// Do updates in parallel.
	for (int iCell = 0; iCell<(Xcells); iCell++){	
		for(int jCell = 0; jCell<(Ycells); jCell++){
		
			if (jCell == 0){ 		// Bottom reflecting wall.
				if (iCell == 0){ 	// Left corner (3 neighbours)
					Tnew[iCell][jCell] = C[iCell][jCell].TOut() + 2./3.*dt*d*(
					C[iCell+1][jCell].TOut() + 
					C[iCell][jCell+1].TOut() +	
					C[iCell+1][jCell+1].TOut() - 
					3.*C[iCell][jCell].TOut()) - dt/taui*C[iCell][jCell].TOut();
					}
				else if (iCell ==(Xcells-1)){ // Right corner (2 neighbours)
					Tnew[iCell][jCell] = C[iCell][jCell].TOut() + 2./3.*dt*d*(
					C[iCell-1][jCell].TOut() + 
					C[iCell][jCell+1].TOut() - 
					2.*C[iCell][jCell].TOut()) - dt/taui*C[iCell][jCell].TOut();
					}
				else{	// Bottom row (4 neighbours)
					Tnew[iCell][jCell] = C[iCell][jCell].TOut() + 2./3.*dt*d*(
					C[iCell+1][jCell].TOut() + 
					C[iCell-1][jCell].TOut() + 
					C[iCell][jCell+1].TOut() +	
					C[iCell+1][jCell+1].TOut()-
					4.*C[iCell][jCell].TOut()) - dt/taui*C[iCell][jCell].TOut();
					}
				}
				
			else if ((jCell>0)&&(jCell<(Ycells-1))){
				if (iCell == 0){ // Left reflecting boundary/wall
					if ((jCell % 2) == 0){ // (5 neighbours)
						Tnew[iCell][jCell] = C[iCell][jCell].TOut() + 2./3.*dt*d*(
						C[iCell+1][jCell].TOut() + 
						C[iCell][jCell+1].TOut() +	
						C[iCell][jCell-1].TOut() +
						// If the rows are even. 
						C[iCell+1][jCell+1].TOut() +
						C[iCell+1][jCell-1].TOut() - 
						5.*C[iCell][jCell].TOut()) - dt/taui*C[iCell][jCell].TOut();
					}
					else{ // (3 neighbours)
						Tnew[iCell][jCell] = C[iCell][jCell].TOut() + 2./3.*dt*d*(
						C[iCell+1][jCell].TOut() + 
						C[iCell][jCell+1].TOut() +	
						C[iCell][jCell-1].TOut() - 
						3.*C[iCell][jCell].TOut()) - dt/taui*C[iCell][jCell].TOut();
					}
				}
				
				
				else if (iCell ==(Xcells-1)){ // Right boundary/wall
					if ((jCell % 2) == 0){ // (3 neighbours)
						Tnew[iCell][jCell] = C[iCell][jCell].TOut() + 2./3.*dt*d*(
						C[iCell-1][jCell].TOut() + 
						C[iCell][jCell+1].TOut() +	
						C[iCell][jCell-1].TOut() -
						3.*C[iCell][jCell].TOut()) - dt/taui*C[iCell][jCell].TOut();
					}
					else{ 	// (5 neighbours)
						Tnew[iCell][jCell] = C[iCell][jCell].TOut() + 2./3.*dt*d*(
						C[iCell-1][jCell].TOut() + 
						C[iCell][jCell+1].TOut() +	
						C[iCell][jCell-1].TOut() +
						// If the rows are uneven.
						C[iCell-1][jCell+1].TOut() +
						C[iCell-1][jCell-1].TOut() - 
						5.*C[iCell][jCell].TOut()) - dt/taui*C[iCell][jCell].TOut();
					}
				}
						
				else{ 		// Interior (6 neighbours)
					if ((jCell % 2) == 0){
						Tnew[iCell][jCell] = C[iCell][jCell].TOut() + 2./3.*dt*d*(
						C[iCell+1][jCell].TOut() + 
						C[iCell-1][jCell].TOut() + 
						C[iCell][jCell+1].TOut() +	
						C[iCell][jCell-1].TOut() +
						// If the rows are even. 
						C[iCell+1][jCell+1].TOut() +
						C[iCell+1][jCell-1].TOut() - 
						6.*C[iCell][jCell].TOut()) - dt/taui*C[iCell][jCell].TOut();
					}
			
					else {
						Tnew[iCell][jCell] = C[iCell][jCell].TOut() + 2./3.*dt*d*(
						C[iCell+1][jCell].TOut() + 
						C[iCell-1][jCell].TOut() + 
						C[iCell][jCell+1].TOut() +	
						C[iCell][jCell-1].TOut() +
						// If the rows are uneven.
						C[iCell-1][jCell+1].TOut() +
						C[iCell-1][jCell-1].TOut() - 
						6.*C[iCell][jCell].TOut()) - dt/taui*C[iCell][jCell].TOut();
					}	
				}
			}
			
			else{ // Upper absorbing wall
				if (iCell == 0){ // (2 neighbours)
					Tnew[iCell][jCell] = C[iCell][jCell].TOut() + 2./3.*dt*d*(
					C[iCell+1][jCell].TOut() +	
					C[iCell][jCell-1].TOut() -
					2.*C[iCell][jCell].TOut()) - dt/taui*C[iCell][jCell].TOut();
				}
				else if (iCell ==(Xcells-1)){ // (3 neighbours)
					Tnew[iCell][jCell] = C[iCell][jCell].TOut() + 2./3.*dt*d*( 
					C[iCell-1][jCell].TOut() + 
					C[iCell][jCell-1].TOut() +
					C[iCell-1][jCell-1].TOut() - 
					3.*C[iCell][jCell].TOut()) - dt/taui*C[iCell][jCell].TOut();
				}
				// Cells are not updated and will act as an open boundary.
				
				/*else{ // (4 neighbours)
					Tnew[iCell][jCell] = C[iCell][jCell].TOut() + 2./3.*dt*d*(
					C[iCell+1][jCell].TOut() + 
					C[iCell-1][jCell].TOut() + 
					C[iCell][jCell-1].TOut() +
					C[iCell-1][jCell-1].TOut() - 
					4.*C[iCell][jCell].TOut())
					- dt/taui*C[iCell][jCell].TOut();
				}*/
			}
			

		// Set negative values of T to 0.
		if (Tnew[iCell][jCell] < 0){Tnew[iCell][jCell]=0;}
		}
	}

	// Reassign the diffused T values
	# pragma omp parallel for private (jCell)
	for (int iCell = 0; iCell<(Xcells); iCell++){
		for	(int jCell = 0; jCell<(Ycells); jCell++){
			C[iCell][jCell].ReassignT(Tnew[iCell][jCell]);}}

	/* Since the even rows are shifted and the distance between all of the cells should be the same, dx,
   the x- and y-distance has to be modified. For the even rows, 0.5 is added to the x-position. The y-distance
   between each row is changed to sqrt(3./4.) instead of 1. This information is contained in the MATLAB-file.*/

	t = t + dt;
	}

	myfile.close();
	myfile1.close();
	myfile2.close();
	myfile3.close();
	myfile4.close();
	myfile5.close();

	ofstream myfile6;
	string fileName6;
	fileName6 = "DATA_time_"+nT+ "_S_"+nS+"_" + ens +".txt";
	myfile6.open((filepath+fileName6).c_str());
	myfile6 << Xcells << " " << Ycells << " " << tstart << " " << tstop << " " << PrintDummy 
		<< " " << newS1 << " ";
	myfile6.close();
	
	double stop = omp_get_wtime();
	printf("Run time: %f\n", stop-start);
	
return 0;
}