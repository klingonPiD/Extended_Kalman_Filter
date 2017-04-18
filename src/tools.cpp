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
	float minThresh = 0.0001;
	//recover state parameters
	float px = (fabs(x_state(0)) < 0.001) ? 0.001 : x_state(0);
	float py = (fabs(x_state(1)) < 0.001) ? 0.001 : x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float p2 = px*px + py*py;
	if (p2 <= minThresh)
	{
		cout << "Warning - Tools::CalculateJacobian - low px/py values" << endl;
		p2 = minThresh;
	}

	float den = sqrt(p2);	
	float dpdpx = px / den;
	float dpdpy = py / den;
	float dpdvx = 0.0;
	float dpdvy = 0.0;

	den = p2;
	float dphidpx = -py / den;
	float dphidpy = px / den;
	float dphidvx = 0.0;
	float dphidvy = 0.0;

	den = pow((pow(px, 2) + pow(py, 2)), 3.0 / 2.0);
	float dprdpx = (vx*py - vy*px)*py / den;
	float dprdpy = (vy*px - vx*py)*px / den;
	float dprdvx = dpdpx;
	float dprdvy = dpdpy;

	//compute the Jacobian matrix
	Hj << dpdpx, dpdpy, dpdvx, dpdvy,
		dphidpx, dphidpy, dphidvx, dphidvy,
		dprdpx, dprdpy, dprdvx, dprdvy;

	return Hj;
}
