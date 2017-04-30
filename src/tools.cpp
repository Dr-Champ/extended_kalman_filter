#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO: Done !
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0.0, 0.0, 0.0, 0.0;

  //accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){
	VectorXd residual = estimations[i] - ground_truth[i];
	
	// coefficient-wise multiplication
	residual = residual.array() * residual.array();
	rmse += residual;
  }

  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO: Done !
    * Calculate a Jacobian here.
  */
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float c1 = px*px+py*py;

  //check division by zero
  if (fabs(c1) < 0.0001) {
  	/*Hj << 0, 0, 0, 0,
  		  0, 0, 0, 0,
  		  0, 0, 0, 0;
  	return Hj;*/
  	c1 = 0.0001;
  }

  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  MatrixXd Hj = MatrixXd(3, 4);
  Hj << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		py * (vx * py - vy * px) / c3, px * (px * vy - py * vx) / c3, px / c2, py / c2;
  return Hj;
}
