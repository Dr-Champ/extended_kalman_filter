#include "kalman_filter.h"
#include <math.h>
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  /**
  TODO: Done !
    * predict the state
  */
  x_ = F_ * x_;
  P_ = (F_ * P_ * F_.transpose()) + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO: Done !
    * update the state by using Kalman Filter equations
  */
  VectorXd y_ = z - (H_ * x_);
  MatrixXd S_ = (H_ * P_ * H_.transpose()) + R_;
  MatrixXd K_ = P_ * H_.transpose() * S_.inverse();
  x_ = x_ + (K_ * y_);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K_ * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO: Done !
    * update the state by using Extended Kalman Filter equations
  */
  double ro_ = sqrt(pow(x_(0), 2) + pow(x_(1), 2));
  double phi_ = atan2(x_(1), x_(0));
  double ro_dot_ = (x_(0) * x_(2) + x_(1) * x_(3)) / ro_;
  VectorXd z_ = VectorXd(3);
  z_ << ro_,
        phi_,
        ro_dot_;
  VectorXd y_ = z - z_;

  // ensure proper phi value
  if (y_(1) < -M_PI) {
    y_(1) = y_(1) + (2 * M_PI);
  } else if (y_(1) > M_PI) {
    y_(1) = y_(1) - (2 * M_PI);
  }

  MatrixXd S_ = (H_ * P_ * H_.transpose()) + R_;
  MatrixXd K_ = P_ * H_.transpose() * S_.inverse();
  x_ = x_ + (K_ * y_);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K_ * H_) * P_;
}










