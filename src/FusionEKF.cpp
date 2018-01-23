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
FusionEKF::FusionEKF() 
{
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  P_       = MatrixXd(4, 4);
  F_       = MatrixXd(4, 4);
  Q_       = MatrixXd(4, 4);
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_  << 1, 0, 0, 0,
	       0, 1, 0, 0;

  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
	0, 0, 10, 0,
	0, 0, 0, 10;

  F_ << 1, 0, 1,  0,
        0, 1, 0,  1,
        0, 0, 1,  0,
        0, 0, 0,  1;

  Q_ << 0, 0, 0,  0,
        0, 0, 0,  0,
        0, 0, 0,  0,
        0, 0, 0,  1;

  noise_ax = 5;
  noise_ay = 5;  
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) 
{
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  
   if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    VectorXd x_l_init(2);
    VectorXd x_r_init(3);
    VectorXd x_init(4);
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
    {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_r_init = measurement_pack.raw_measurements_;
      x_init << x_r_init(0)*cos(x_r_init(1)), x_r_init(0)*sin(x_r_init(1)), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      /**
      Initialize state.
      */
      x_l_init = measurement_pack.raw_measurements_;
      x_init << x_l_init(0), x_l_init(1), 0, 0;
    }
    ekf_.Init(x_init, P_, F_, H_laser_, R_laser_, R_radar_, Q_);
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
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
   float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
   previous_timestamp_ = measurement_pack.timestamp_;
   
   float dt_2 = dt * dt;
   float dt_3 = dt_2 * dt;
   float dt_4 = dt_3 * dt;
 
   //update state transition matrix F
   ekf_.F_ << 1, 0, dt, 0,
              0, 1, 0,  dt,
              0, 0, 1,  0,
              0, 0, 0,  1;

   //update the process covariance matrix Q
   ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
	       0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
	       dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
	       0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

   ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
  {
    // Radar updates
    // Calculate Hj based on produced states
    Hj_ = tools.CalculateJacobian(ekf_.x_, Hj_);
    ekf_.Update(measurement_pack.raw_measurements_, Hj_);
  } 
  else 
  {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }
}
