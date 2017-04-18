#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	unsigned long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	//Get jacobian Hj
	double minThresh = 0.0001;
	Tools tools;
	MatrixXd Hj = tools.CalculateJacobian(x_);

	//compute h_x - in a seperate func
	VectorXd h_x(3);
	h_x << 0.0, 0.0, 0.0;
	//recover state parameters
	double px = (fabs(x_(0)) < minThresh) ? minThresh : x_(0);
	double py = (fabs(x_(1)) < minThresh) ? minThresh : x_(1);
	double vx = x_(2);
	double vy = x_(3);

	//compute h_x
	h_x[0] = sqrt(px*px + py*py);
	if(h_x[0] < minThresh)
	{
		cerr << "Error - KalmanFilter::UpdateEKF - division by zeros" << endl;
	}
	h_x[1] = atan2(py, px);
	h_x[2] = (px*vx + py*vy) / h_x[0];


	VectorXd z_pred = h_x;// H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Hjt = Hj.transpose();
	MatrixXd S = Hj * P_ * Hjt + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHjt = P_ * Hjt;
	MatrixXd K = PHjt * Si;

	//new estimate
	//Conversion back to polar co-ordinates not required since Kalman gain K
	//implicitly handles this for us
	x_ = x_ + (K * y);
	unsigned long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj) * P_;

}
