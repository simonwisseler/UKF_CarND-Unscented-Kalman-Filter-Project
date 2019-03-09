#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;


/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // state dimension
    n_x_ = 5;

    // initial state vector
    x_ = VectorXd(n_x_);

    // initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);
    P_ <<   0.1, 0, 0, 0, 0,
            0, 0.1, 0, 0, 0,
            0, 0, 0.2, 0, 0,
            0, 0, 0, 1.0, 0,
            0, 0, 0, 0, 1.0;

    // process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 1.5;

    // process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.3;

    // augmented state dimension
    n_aug_ = n_x_ + 2;
    
    // sigma point spreading parameter --> heuristic from Julier, S., and J. Uhlmann. 1997. A new extension of the Kalman filter to nonlinear systems
    lambda_ = 3 - n_x_;
    
    // weights of sigma points
    weights_ = VectorXd(2*n_aug_ + 1);
    weights_.fill(0.5/(n_aug_ + lambda_));
    weights_(0) = lambda_/(lambda_ + n_aug_);

    // predicted sigma points matrix
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    
    // laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    // radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;

    // laser measurement noise covariance matrix
    R_laser_ = MatrixXd(2, 2);
    R_laser_ << std_laspx_*std_laspx_,0,
                0,std_laspy_*std_laspy_;
    
    // radar measurement noise covariance matrix
    R_radar_ = MatrixXd(3, 3);
    R_radar_ << std_radr_*std_radr_, 0, 0,
                0, std_radphi_*std_radphi_, 0,
                0, 0,std_radrd_*std_radrd_;
    
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    
    if (!is_initialized_) {
        
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
        } else {
            double rho = meas_package.raw_measurements_[0]; // range
            double phi = meas_package.raw_measurements_[1]; // bearing
            double rho_dot = meas_package.raw_measurements_[2]; // range rate (https://en.wikipedia.org/wiki/Range_rate)
            
            double x = rho * cos(phi);
            double y = rho * sin(phi);
            // best estimate for v as component of velocity normal to rho not captured by radar measurement
            double v = rho_dot
            x_ << x, y, v, 0, 0;
        }
        previous_timestamp_ = meas_package.timestamp_ ;
        is_initialized_ = true;
        return;
    }
    
    float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; // time (in sec) elapsed between current and previous measurement
    previous_timestamp_ = meas_package.timestamp_;
    
    Prediction(dt);
    
    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
        UpdateLidar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && user_radar_) {
        UpdateRadar(meas_package);
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}
