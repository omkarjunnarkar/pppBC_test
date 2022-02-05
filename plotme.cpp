#include<iostream>
#include<vector>
#include"src/rapidcsv.h"
#include"gnuplot-iostream.h"
#include"plotme.h"

using namespace std;
using namespace rapidcsv;

int plot() {

	Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");
	Document dy("displacement_2_bar.csv");
	Document dstress("stress.csv");
	Document dstrain("strain.csv");
	vector<double> getdisplacement = dy.GetColumn<double>(0);
	vector<double> getstress = dstress.GetColumn<double>(0);
	vector<double> getstrain = dstrain.GetColumn<double>(0);

	int data_size = getstress.size();
	//cout << "data size= " << data_size << endl;;
	vector<vector<double>>stress_strain(data_size,vector<double>(2));

	for (int i = 0; i < data_size; i++) {

		stress_strain[i][0] = getstrain[i];
		stress_strain[i][1] = getstress[i];
		//cout << getstress[i] << " , " << getstrain[i] << endl;	
		
	}

	gp << "set multiplot layout 1,2 rowsfirst\n";
	gp << "set xlabel 'Strain' font 'Times - Roman, 10'\n";
	gp << "set ylabel 'Stress' font 'Times - Roman, 10'\n";
	gp << "set xtics font 'Arial, 7'\n";
	gp << "set ytics font 'Arial, 7'\n";

	gp << "set title 'Stress-Strain Curve'\n";
	gp << "plot '-' with lines title 'Stress'\n";

	gp.send1d(stress_strain);

	//gp << "set size 0.36, 0.5\n";
	gp << "set title 'Displacement Curve'\n";
	gp << "set xlabel 'Time-Steps' font 'Times - Roman, 10'\n";
	gp << "set ylabel 'Displacement' font 'Times - Roman, 10'\n";
	gp << "set xtics font 'Arial, 7'\n";
	gp << "set ytics font 'Arial, 7'\n";

	gp << "plot '-' with lines title 'u'\n";

	gp.send1d(getdisplacement);

	cin.get();

};