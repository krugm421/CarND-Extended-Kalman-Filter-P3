#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  //measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  //measurement matrix lidar
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //jacobian radar
  Hj_ = MatrixXd(3, 4);
  Hj_ << 0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  //create and initialize kalman filter 
  //ekf_ = KalmanFilter();

  //VectorXd x_ = VectorXd(4);
  //x_ << 0,0,0,0;

  //MatrixXd F_ = MatrixXd(4,4);
  //F_ << 1, 0, 1, 0,
  //      0, 1, 0, 1,
  //      0, 0, 1, 0,
  //      0, 0, 0, 1;

  //MatrixXd Q_ = MatrixXd(4, 4);
  //Q_ << 0, 0, 0, 0,
  //      0, 0, 0, 0,
  //      0, 0, 0, 0,
  //      0, 0, 0, 0;

  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ <<  1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;

  // MatrixXd H_ = MatrixXd(2, 4); 
  // H_ << 1, 0, 0, 0,
  //       0, 1, 0, 0;

  // MatrixXd R_ = MatrixXd(2,2);
  // R_ << 0.0225, 0,
  //       0, 0.0225;

  // ekf_.Init(x_, P_, F_, H_, R_, Q_);




}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0, 0, 0, 0;

    if(measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    else if(measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      //Convert radar from polar to cartesian coordinates 
      float ro = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];
      float ro_dot = measurement_pack.raw_measurements_[2];

      float x = ro * cos(theta);
      x = (x < 0.0001) ? 0.0001 : x;
      float y = ro * sin(theta);
      y = (y < 0.0001) ? 0.0001 : y;
      float v_x = ro_dot * cos(theta);
      float v_y = ro_dot * sin(theta);

      ekf_.x_ << x, y, v_x, v_y;
      //ekf_.x_ << ro * sin(theta), ro * cos(theta), ro_dot * sin(theta), ro_dot * cos(theta);  // <-- initialize without the speed measurement from the radar measurement
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  //update the transition matrix
  float delta_t = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0; 
  previous_timestamp_ = measurement_pack.timestamp_;
  
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, delta_t, 0,
             0, 1, 0, delta_t,
             0, 0, 1, 0,
             0, 0, 0, 1;

  //update the process noise covariance matrix Q
  double noise_ax = 9.0;
  double noise_ay = 9.0;

  double dt_2 = delta_t * delta_t; //dt^2
  double dt_3 = dt_2 * delta_t; //dt^3
  double dt_4 = dt_3 * delta_t; //dt^4
  double dt_4_4 = dt_4 / 4; //dt^4/4
  double dt_3_2 = dt_3 / 2; //dt^3/2
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4_4 * noise_ax, 0, dt_3_2 * noise_ax, 0,
	         0, dt_4_4 * noise_ay, 0, dt_3_2 * noise_ay,
	         dt_3_2 * noise_ax, 0, dt_2 * noise_ax, 0,
 	         0, dt_3_2 * noise_ay, 0, dt_2 * noise_ay;
  

  ekf_.Predict();

  /**
  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
