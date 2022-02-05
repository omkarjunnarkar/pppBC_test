#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<fstream>
#include<vector>
#include"src/rapidcsv.h"
#include"material.h"
#include"element.h"
#include"plotme.h"

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

using namespace std;
using namespace Eigen;
using namespace rapidcsv;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

void main() {

	cout << "Starting the Finite Element Program" << endl;
	//Creating csv files for storing stress, train and displacement data

	ofstream mydisplacefile("displacement_2_bar.csv");
	ofstream mystressfile("stress.csv");
	ofstream mystrainfile("strain.csv");
	ofstream Allstrainfile("All_strain.csv");
	ofstream Allstressfile("All_stress.csv");
	ofstream Residualfile("Residual.csv");
	ofstream input("InputFileFEA2Bar.txt");
	ofstream bmat("B-Matrix.txt");
	ofstream amat("Assignment-Matrix.txt");
	ofstream ar_len("Areas-Lengths-After-Meshing.txt");
	ofstream loadhist("LoadingSteps.csv");
	ofstream ReactionForces("ReactionForces.csv");
	

	// Reading User Input File

	Document inp("UserInputFile2BarFEA.csv");
	vector<double> user_input = inp.GetColumn<double>(0);

	double ForceMax = user_input[0];			//Maximum Force
	double Area_1 = user_input[1];					//csArea of left element
	double Area_2 = user_input[2];				//csArea of right element
	double Length_1 = user_input[3];				//Length of left element
	double Length_2 = user_input[4];				//Length of right element
	double E = user_input[5]*1e3;				//Youngs Modulus
	int steps_num = 10000;			//Number of Steps to solve the problem
	double steps = steps_num;
	double t_tot = user_input[6];				//Total time for which force is applied
	double yield = user_input[7];			//Yield Strength
	double eta = user_input[8];					//Viscosity
	int elements = user_input[9];				//Number of elements in which problem is to be divided

	//double ForceMax = 4700;			//Maximum Force
	//double Area_1 = 6;					//csArea of left element
	//double Area_2 = 12;				//csArea of right element
	//double Length_1 = 30;				//Length of left element
	//double Length_2 = 60;				//Length of right element
	//double E = 120000;				//Youngs Modulus
	//int steps_num = 10000;			//Number of Steps to solve the problem
	//double steps = steps_num;
	//double t_tot = 0.025;				//Total time for which force is applied
	//double yield = 240;			//Yield Strength
	//double eta = 20;					//Viscosity
	//int elements = 2;				//Number of elements in which problem is to be divided
	double force = 0;					
	int gauss_weights = 2;

	MatrixXd Areas, Lengths;

	input << "Input Data for FEA Code: " << endl;
	input << "Max. Force (Newton) = " << ForceMax << endl;
	input << "Area of Bar 1 (mm^2) = " << Area_1 << endl;
	input << "Area of Bar 2 (mm^2) = " << Area_2 << endl;
	input << "Length of Bar 1 (mm) = " << Length_1 << endl;
	input << "Length of Bar 1 (mm) = " << Length_2 << endl;
	input << "Young's Modulus (GPa) = " << E/1e3 << endl;
	input << "Yield Stress (MPa) = " << yield << endl;
	input << "Viscosity (1/s) = " << eta << endl;
	input << "Total Time (s) = " << t_tot << endl;
	input << "Time Steps = " << steps << endl;
	input << "No. of Elements = " << elements << endl;
	input << "Gauss Points, Weights = " << "1 , " << gauss_weights << endl;

	steps = t_tot / double(steps);

	//Calling Element Routine to get Areas & Lengths of all the elements

	tie(Areas, Lengths) = get_areas_lengths(elements, Area_1, Area_2, Length_1, Length_2);	
	
	input << "\nList of Areas after Seeding: \n" << Areas << endl;
	input << "\nList of Lengths after Seeding: \n" << Lengths << endl;

	MatrixXd F_ext = MatrixXd::Zero(elements + 1,1);
	MatrixXd E_mat= E*MatrixXd::Ones(elements , 1);
	MatrixXd eps_pl = MatrixXd::Zero(elements , 1);
	MatrixXd epsilon = MatrixXd::Zero(elements , 1);
	MatrixXd u = MatrixXd::Zero(elements + 1, 1);
	MatrixXd Kt, F_int, strain, stress;
	MatrixXd temp_mat1,temp_mat2;
	
	int count = 0;
	
	for (int e = 0; e < elements; e++) {
		bmat<<"B - Matrix for Element ["<<e<<"] of length = " << Lengths(e,0) << " mm is : \n" << B_mat(Lengths(e, 0)) << endl;
		amat << "Assignment Matrix for Element [" << e << "] in system of " << elements << " Elements is :\n" << assign(e, elements) << endl;
		ar_len << "Area of Element [" << e << "] = " << Areas(e, 0) << " Length of Element [" << e << "] = " << Lengths(e, 0) << endl;
	}
	bmat.close();
	amat.close();
	ar_len.close();

	//Iterating through all force steps

	while ((abs(force) + (steps / t_tot) * abs(ForceMax))<= abs(ForceMax) + (steps / t_tot) * abs(ForceMax)) {
		if (count != 0) { force = force + (steps / t_tot) * ForceMax; }		//Force Increment
		
		count++;
		
		//cout << "Force = " << force << endl;

		/*
		
		For ORIGINAL & CASE-1 :
		F_ext(int(elements / 2), 0) = force;	//Assigning the node to apply force in Global External Force Matrix
		
		*/

		/*

		For CASE-2 & CASE-3:
		F_ext(elements , 0) = force;	//Assigning the node to apply force in Global External Force Matrix

		*/

		/*

		For CASE-4 :
		F_ext(elements , 0) = -force;	//Assigning the node to apply force in Global External Force Matrix

		*/

		/*

		For CASE-5 :
		F_ext(int(elements/2) , 0) = -force;	//Assigning the node to apply force in Global External Force Matrix

		*/

		F_ext(int(elements / 2), 0) = -force;	//Assigning the node to apply force in Global External Force Matrix

		loadhist << force << endl;
		//Entering the Newton-Raphson Method

		for (int NR = 0; NR < 6; NR++) {

			//Calling Material Routine to get the Global Tangent Stiffness Matrix, Global Internal Force Matrix, Plastic Strain, Total Strain & Stress Values of All Elements

			tie(Kt, F_int, eps_pl, strain, stress) = material_routine(elements, gauss_weights, Lengths, Areas, E_mat, u, steps, eta, eps_pl, yield, E);
			
			//cout << "\nKt_Global = \n" << Kt << endl;
			
			/*
			
			For BC: Both Ends Fixed, Force on Center Node : ORIGINAL
			
			MatrixXd Kt_red = Kt.block(1, 1, Kt.cols() - 2, Kt.rows() - 2);
			MatrixXd G_red(elements - 1, 1),u_red(elements-1,1);					//Reducing the Residual
			
			for (int c = 1; c < elements; c++) { G_red(c - 1, 0) = G(c, 0); }
			for (int c = 1; c < elements; c++) { u_red(c - 1, 0) = u(c, 0); }

			for (int c = 1; c < elements; c++) { u(c, 0) = u_red(c - 1, 0); }

			*/

			/*
			
			For BC: Second end Free, Force on Center Node : CASE-1 & CASE-2 & CASE-3 & CASE-4 & CASE-5

			MatrixXd Kt_red = Kt.block(1, 1, Kt.cols() - 1, Kt.rows() - 1);
			MatrixXd G_red(elements , 1),u_red(elements,1);						//CHANGE FOR BC		//Reducing the Residual

			for (int c = 1; c < elements+1; c++) { G_red(c - 1, 0) = G(c, 0); }	//CHANGE FOR BC
			for (int c = 1; c < elements+1; c++) { u_red(c - 1, 0) = u(c, 0); }	//CHANGE FOR BC

			for (int c = 1; c < elements+1; c++) { u(c, 0) = u_red(c - 1, 0); }	//CHANGE FOR BC

			*/

			MatrixXd Kt_red = Kt.block(1, 1, Kt.cols() - 1, Kt.rows() - 1);		//CHANGE FOR BC		//Reducing the Tangent Stiffness
			MatrixXd inv_Kt_red = Kt_red.inverse();

			MatrixXd G = F_int - F_ext;

			//cout << "G=" << G << endl;

			MatrixXd G_red(elements , 1),u_red(elements,1);						//CHANGE FOR BC		//Reducing the Residual 
			for (int c = 1; c < elements+1; c++) { G_red(c - 1, 0) = G(c, 0); }	//CHANGE FOR BC
			for (int c = 1; c < elements+1; c++) { u_red(c - 1, 0) = u(c, 0); }	//CHANGE FOR BC

			MatrixXd del_u_red = inv_Kt_red * G_red;								//Computing Delta U
			
			u_red = u_red - del_u_red;
			for (int c = 1; c < elements+1; c++) { u(c, 0) = u_red(c - 1, 0); }	//CHANGE FOR BC
			
			if (NR == 5) { 
				
				Residualfile << G_red.norm() << endl;
				ReactionForces << F_int(0, 0) << " , " << F_int(elements , 0) << endl;
				/*if (del_u_red.norm() < 0.005 * u_red.norm()) {
					cout << "NR Converged" << endl;
				}
				else cout << "NR NOT Converged" << endl;*/
			}
			
		}
		
		//Appending values of Center Node/First Element to the csv files

		mydisplacefile << u((elements / 2) , 0) << endl;
		mystressfile << stress((elements / 2)-1, 0) << endl;
		mystrainfile << strain((elements / 2) - 1, 0) << endl;

		temp_mat1 = strain.reshaped(1, elements);
		temp_mat2 = stress.reshaped(1, elements);
		for (int fc = 0; fc < elements; fc++) {
			Allstrainfile << temp_mat1(0,fc) << "," ;
			Allstressfile << temp_mat2(0, fc) << "," ;
		}
		Allstrainfile << endl;
		Allstressfile << endl;
		
	}

	mydisplacefile.close();
	mystressfile.close();
	mystrainfile.close();
	Allstrainfile.close();
	Allstressfile.close();
	Residualfile.close();
	input.close();
	loadhist.close();
	ReactionForces.close();
	
	cout << "Computation Done ! " << endl;
	cout << "Activating GNUPLOT using vcpkg and gnupuplot-iostream..\n";
	plot();

}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
