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

  /**
  TODO: Done !
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  ekf_.x_ = VectorXd(4);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd(4, 4);

  // Initialize P_ such that it represents my uncertainty about locations, and
  // my even more uncertain about the velocities.
  ekf_.P_ <<  100, 0, 0, 0,
              0, 100, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;
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
    /**
    TODO: Done !
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    //ekf_.x_ = VectorXd(4);
    //ekf_.x_ << 1, 1, 1, 1;
    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);

      // only initialize location information, not velocity. 
      ekf_.x_ << (ro * cos(phi)), (ro * sin(phi)), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      float x = measurement_pack.raw_measurements_(0);
      float y = measurement_pack.raw_measurements_(1);

      // only initialize location information, not velocity.
      ekf_.x_ << x, y, 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO: Done !
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  double delta_t = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; // s
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_ <<  1, 0, delta_t, 0,
              0, 1, 0, delta_t,
              0, 0, 1, 0,
              0, 0, 0, 1;

  MatrixXd G = MatrixXd(4, 2);
  double g1 = delta_t * delta_t / 2.0;
  double g2 = delta_t;
  G <<  g1, 0,
        0, g1,
        g2, 0,
        0, g2;
  MatrixXd Q_a = MatrixXd(2, 2);
  Q_a <<  9, 0,
          0, 9;
  ekf_.Q_ = G * Q_a * G.transpose();
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO: Done !
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    float ro = measurement_pack.raw_measurements_(0);
    float phi = measurement_pack.raw_measurements_(1);
    float ro_dot = measurement_pack.raw_measurements_(2);
    VectorXd z = VectorXd(3);
    z << ro,
         phi,
         ro_dot;
    Hj_ = tools.CalculateJacobian(ekf_.x_); // Jacobian, evaluated at predicted x
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(z);
  } else {
    // Laser updates
    float x = measurement_pack.raw_measurements_(0);
    float y = measurement_pack.raw_measurements_(1);
    VectorXd z = VectorXd(2);
    z << x,
         y;
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}










