#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<fstream>
#include"element.h"

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

using namespace std;
using namespace Eigen;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

MatrixXd B_mat(double le) {
	MatrixXd b(1, 2);
	b(0, 0) = -1 / le;
	b(0, 1) = 1 / le;
	//cout << "B-Matrix for Element of length = " << le << " is : \n" << b << endl;
	return b;
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

int signum(double x) {
	if (x >= 0) return 1;
	else return -1;
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*To compute the assignment matrix for each element*/
MatrixXd assign(int ele_index_cpp, int numele) {

	//ele_index is current element number, not defined here!! different for each ele !!

	MatrixXd a_element = MatrixXd::Zero(2, numele + 1);
	a_element(0, ele_index_cpp) = 1;
	a_element(1, ele_index_cpp + 1) = 1;
	return a_element;
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

tuple<MatrixXd, MatrixXd>get_areas_lengths(int elements, double ar_1, double ar_2, double len_1, double len_2) {
	MatrixXd areas(elements, 1);
	MatrixXd lengths(elements, 1);
	for (int i = 0; i < elements; i++) {
		if (i < elements / 2) {
			areas(i, 0) = ar_1;
			lengths(i, 0) = len_1 / (elements / 2);
		}
		else {
			areas(i, 0) = ar_2;
			lengths(i, 0) = len_2 / (elements - (elements / 2));
		}
	}
	
	return { areas,lengths };
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/



