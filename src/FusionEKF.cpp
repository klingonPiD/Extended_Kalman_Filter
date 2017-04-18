#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
	  0, 1, 0, 0;

  // process noise
  P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 1, 0,
	  0, 0, 0, 1;
  //measurement noise
  noise_ax_ = 9;
  noise_ay_ = 9;
 
  //rest of the stuff
  F_ = MatrixXd(4, 4);
  Q_ = MatrixXd(4, 4);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
		x_ = VectorXd(4);
		VectorXd z_meas = VectorXd(3);
		z_meas << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), measurement_pack.raw_measurements_(2);
		double pho = z_meas(0);
		double phi = z_meas(1);
		double px = pho * cos(phi);
		double py = pho * sin(phi);
		x_ <<px, py, 0.0, 0.0;
		ekf_.Init(x_, P_, F_, Hj_, R_radar_, Q_);
		previous_timestamp_ = measurement_pack.timestamp_;
		is_initialized_ = true;
		return;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
		x_ = VectorXd(4);
		x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0.0, 0.0;
		ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);
		previous_timestamp_ = measurement_pack.timestamp_;
		is_initialized_ = true;
		return;
    }

    // done initializing, no need to predict or update
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  //compute the time elapsed between the current and previous measurements
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  // TODO: YOUR CODE HERE
  //1. Modify the F matrix so that the time is integrated
  ekf_.F_ << 1, 0, dt, 0,
	  0, 1, 0, dt,
	  0, 0, 1, 0,
	  0, 0, 0, 1;
  //2. Update the process noise covariance matrix
  //2. Set the process covariance matrix Q
  double noise_ax = noise_ax_;// pow(noise_ax_, 2);
  double noise_ay = noise_ay_;// pow(noise_ay_, 2);
  double dt4 = pow(dt, 4);
  double dt3 = pow(dt, 3);
  double dt2 = pow(dt, 2);
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << (dt4 / 4.0)*noise_ax, 0, (dt3 / 2.0)*noise_ax, 0,
	  0, (dt4 / 4.0)*noise_ay, 0, (dt3 / 2.0)*noise_ay,
	  (dt3 / 2.0)*noise_ax, 0, dt2*noise_ax, 0,
	  0, (dt3 / 2.0)*noise_ay, 0, dt2*noise_ay;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  ekf_.H_ = Hj_;
	  ekf_.R_ = R_radar_;
	  VectorXd z_meas(3);
	  z_meas << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), measurement_pack.raw_measurements_(2);
	  ekf_.UpdateEKF(z_meas);
  } else {
    // Laser updates
	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  VectorXd z_meas(2);
	  z_meas << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1);
	  ekf_.Update(z_meas);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
