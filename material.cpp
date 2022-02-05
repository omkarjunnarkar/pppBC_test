#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<fstream>
#include"material.h"
#include"element.h"

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

using namespace std;
using namespace Eigen;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

double Lambda(double sigma, double yield) {
	if (abs(sigma) <= yield) { return 0; }
	else return ((abs(sigma) / yield) - 1);
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

tuple<double, double, double>get_Stress_Strains(double eps_p_m, double epsilon, double eta, double h, double t_total, double yield, double E) {
	
	int s = (t_total / h);
	double Ct = E;
	double sigma_m,eps_p_next, eps_p_prev, sigma_next;
	
	// 1 step Euler Explicit

	for (int i = 0; i<s; i++) {
		sigma_m = E * (epsilon - eps_p_m);
		eps_p_next = eps_p_m + h * eta * signum(sigma_m) * Lambda(sigma_m, yield);

		// 1000 Steps Euler Implicit

		for (int j = 0; j < 1000; j++) {
			eps_p_prev = eps_p_next;
			sigma_next = E * (epsilon - eps_p_next);
			eps_p_next = eps_p_m + h * eta * signum(sigma_next) * Lambda(sigma_next, yield);

			if (abs(eps_p_next - eps_p_prev) < 1e-5) {
				break;
			}
		}
		
		sigma_m = E * (epsilon - eps_p_next);						//Updating final values of STress and Plastic Strain
		eps_p_m = eps_p_next;

	}

	if (abs(sigma_m) > yield) {
		Ct = (E * yield) / (yield + E * eta * t_total);		//Updating Tangent Stiffness for Plastic Case
		//cout << "Ct = " << Ct << endl;
	}
	else Ct = E;

	//cout << "Tangent Stiffness = " << Ct << endl;

	return { sigma_m, eps_p_next, Ct };//,sigmas, ep_p_vals
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd>material_routine(int elements, double gauss_weights, MatrixXd Lengths, MatrixXd Areas, MatrixXd  E_mat, MatrixXd u, double time_step, double eta, MatrixXd eps_pl, double yield, double E) {

	MatrixXd F_int= MatrixXd::Zero(elements + 1, 1);
	MatrixXd stress_mat = MatrixXd::Zero(elements , 1);
	MatrixXd Kt = MatrixXd::Zero(elements+1 , elements+1);
	MatrixXd epsilon = MatrixXd::Zero(elements , 1);
	MatrixXd u_temp(2, 1);

	for (int ele = 0; ele < elements; ele++) {
		

		u_temp(0, 0) = u(ele,0);
		u_temp(1, 0) = u(ele + 1, 0);
		MatrixXd t = B_mat(Lengths(ele,0)) * u_temp;
		epsilon(ele, 0) = t(0,0);

		double smt, ept, et;
		tie(smt, ept, et) = get_Stress_Strains(eps_pl(ele, 0), epsilon(ele, 0), eta, time_step / 100, time_step, yield, E);
		stress_mat(ele, 0) = smt;
		eps_pl(ele, 0) = ept;
		E_mat(ele, 0) = et;
		
		MatrixXd F_int_temp = B_mat(Lengths(ele, 0)).transpose() * Areas(ele, 0) * stress_mat(ele, 0) * (Lengths(ele, 0) / 2) * gauss_weights;
		MatrixXd Kt_temp = B_mat(Lengths(ele, 0)).transpose() * B_mat(Lengths(ele, 0)) * E_mat(ele, 0) * Areas(ele, 0) * (Lengths(ele, 0) / 2) * gauss_weights;
		/*if (E_mat(0, 0) != E) { 
			cout << "Kt_Element[" << ele << "] = \n" << Kt_temp << endl; }*/
		
		//Converting from Local to Global at one Gauss Point

		F_int(ele, 0) = F_int(ele, 0) + F_int_temp(0, 0);
		F_int(ele + 1, 0) = F_int(ele + 1, 0) + F_int_temp(1, 0);

		MatrixXd A = assign(ele, elements);
		Kt = Kt + (A.transpose() * Kt_temp * A);

	}
	/*if (E_mat(0, 0) != E) { 
		cout << "Kt_Global = \n" << Kt << endl; }*/
	
	return { Kt, F_int, eps_pl, epsilon, stress_mat };

}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/