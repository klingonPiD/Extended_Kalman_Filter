#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if (estimations.size() != ground_truth.size() || estimations.size() < 1)
	{
		cout << "Error - Invalid estimation / ground truth array size";
		return rmse;
	}

	//accumulate error
	for (unsigned int i = 0; i < estimations.size(); i++) {
		VectorXd err = estimations[i] - ground_truth[i];
		err = err.array()*err.array();
		rmse += err;
	}

	//calculate the mean
	rmse = rmse / estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3, 4);
	double minThresh = 0.0001;
	//recover state parameters
	double px = (fabs(x_state(0)) < 0.001) ? 0.001 : x_state(0);
	double py = (fabs(x_state(1)) < 0.001) ? 0.001 : x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	double p2 = px*px + py*py;
	if (p2 <= minThresh)
	{
		cout << "Warning - Tools::CalculateJacobian - low px/py values" << endl;
		p2 = minThresh;
	}

	double den = sqrt(p2);
	double dpdpx = px / den;
	double dpdpy = py / den;
	double dpdvx = 0.0;
	double dpdvy = 0.0;

	den = p2;
	double dphidpx = -py / den;
	double dphidpy = px / den;
	double dphidvx = 0.0;
	double dphidvy = 0.0;

	den = pow(den, 3.0 / 2.0);
	double dprdpx = (vx*py - vy*px)*py / den;
	double dprdpy = (vy*px - vx*py)*px / den;
	double dprdvx = dpdpx;
	double dprdvy = dpdpy;

	//compute the Jacobian matrix
	Hj << dpdpx, dpdpy, dpdvx, dpdvy,
		dphidpx, dphidpy, dphidvx, dphidvy,
		dprdpx, dprdpy, dprdvx, dprdvy;

	return Hj;
}
