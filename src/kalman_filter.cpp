#include "kalman_filter.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &Rl_in, MatrixXd &Rr_in,
                        MatrixXd &Q_in) 
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  Rl_ = Rl_in;
  Rr_ = Rr_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() 
{
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) 
{
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  MatrixXd y;
  MatrixXd S;
  MatrixXd K;
  y = z - H_ * x_;
  S = H_ * P_ * H_.transpose() + Rl_;
  K = P_ * H_.transpose() * S.inverse();

  x_ = x_ + K * y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::Update(const VectorXd &z, const Eigen::MatrixXd &Hj) 
{
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  VectorXd h(3);
  float c1 = sqrt(x_(0)*x_(0) + x_(1)*x_(1));

  h(0) = c1;
  h(1) = atan2(x_(1),x_(0));
  h(2) = (x_(0)*x_(2) + x_(1)*x_(3))/c1;
  if(fabs(c1) < 0.00001)
  {
    h(2) = 0;
  }

  MatrixXd y  = z - h;

  //
  // make sure  -pi < y(1) < pi
  // 
  while (y(1) < -M_PI) 
  { 
      y(1) += 2*M_PI;
  }
  while (y(1) > M_PI) 
  { 
      y(1) -= 2*M_PI;
  }

  MatrixXd S  = Hj * P_ * Hj.transpose() + Rr_;
  MatrixXd K  = P_ * Hj.transpose() * S.inverse();

  x_ = x_ + K * y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;
}
