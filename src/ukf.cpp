/**
 * Unscented Kalman Filter as described in
 * Julier, S., and J. Uhlmann. 1997. A new extension of the Kalman filter to nonlinear systems
 */

#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;


/**
 * Constructor
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
    
    // sigma point spreading parameter --> using value 3 is heuristic from Julier and Uhlmann (1997)
    lambda_ = 3 - n_aug_;
    
    // weights of augmented sigma points
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
    R_laser_ << std_laspx_ * std_laspx_, 0,
                0, std_laspy_ * std_laspy_;
    
    // radar measurement noise covariance matrix
    R_radar_ = MatrixXd(3, 3);
    R_radar_ << std_radr_ * std_radr_, 0, 0,
                0, std_radphi_ * std_radphi_, 0,
                0, 0, std_radrd_ * std_radrd_;
}


/**
 * Destructor
 */
UKF::~UKF() {}


/**
 * ProcessMeasurement Run prediction and update step for radar and/or laser measurements
 * @param meas_package The latest measurement data of either radar or laser
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    
    if (!is_initialized_) {
        
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
        } else {
            double rho = meas_package.raw_measurements_[0]; // range
            double phi = meas_package.raw_measurements_[1]; // bearing
            double rho_dot = meas_package.raw_measurements_[2]; // range rate (https://en.wikipedia.org/wiki/Range_rate)
            
            // transform polar coordinates rho, phi, rho_dot into cartesian coordinates x, y, v
            double x = rho * cos(phi);
            double y = rho * sin(phi);
            
            // best estimate for v as component of velocity normal to rho not captured by radar measurement
            double v = rho_dot;
            
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
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
        UpdateRadar(meas_package);
    }
}


/**
 * Prediction Predicts augmented sigma points, the state, and the state covariance matrix
 * @param delta_t Time between k and k+1 in s
 */
void UKF::Prediction(double delta_t) {
    //generate augmented sigma points
    MatrixXd Xsig_aug = GenerateAugmentedSigmaPoints();

    // predict augmented sigma points
    PredictAugmentedSigmaPoints(Xsig_aug, delta_t);
    
    // predict state mean
    x_.fill(0.0);
    for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }
    
    // predict state covariance
    P_.fill(0.0);
    for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd x_res = Xsig_pred_.col(i) - x_;
        
        NormalizeAngle(x_res, 3); // normalize phi if it is > |2*pi|
        
        P_ = P_ + weights_(i) * x_res * x_res.transpose();
    }
}


/**
 * UpdateLidar Updates the state and the state covariance matrix using a laser measurement
 * @param meas_package The measurement at k+1
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    VectorXd z = meas_package.raw_measurements_;
    
    // dimension of measurement vector (measurement only consists of x and y coordinate)
    int n_z = 2;
    
    // generate sigma points for laser measurement
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig.row(0) = Xsig_pred_.row(0); // px
    Zsig.row(1) = Xsig_pred_.row(1); // py
    
    // predict mean laser measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    
    for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    
    // calculate covariance matrix of predicted laser measurement
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    
    for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    
    // add measurement noise covariance
    S = S + R_laser_;
    
    // calculate cross-correlation
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        //
        VectorXd z_res = Zsig.col(i) - z_pred;
        
        //
        VectorXd x_res = Xsig_pred_.col(i) - x_;
        
        Tc = Tc + weights_(i) * x_res * z_res.transpose();
    }
    
    // calculate Kalman gain
    MatrixXd K = Tc * S.inverse();
    
    // calculate innovation
    VectorXd innovation = z - z_pred;
    
    // update state mean and covariance
    x_ = x_ + K * innovation;
    P_ = P_ - K * S * K.transpose();
    
    // determine normalized innovation squared for testing parameters' goodness of fit
    NIS_laser_ = innovation.transpose() * S.inverse() * innovation;
}


/**
 * UpdateRadar Updates the state and the state covariance matrix using a radar measurement
 * @param meas_package The measurement at k+1
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    VectorXd z = meas_package.raw_measurements_;
    
    // dimension of measurement vector (measurement consists of rho, phi and rho_dot coordinate)
    int n_z = 3;
    
    // generate sigma points for laser measurement
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    
    for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {
        
        double px = Xsig_pred_(0, i);
        double py = Xsig_pred_(1, i);
        double v   = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);
        
        double vx = cos(yaw) * v;
        double vy = sin(yaw) * v;
        
        // measurement model
        Zsig(0, i) = sqrt(px * px + py * py);
        Zsig(1, i) = atan2(py, px);
        Zsig(2, i) = (px * vx + py * vy) / sqrt(px * px + py * py);
    }
    
    // predict mean radar measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    
    for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    
    // calculate covariance matrix of predicted radar measurement
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    
    for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        NormalizeAngle(z_diff, 1);
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    
    // add measurement covariance matrix S
    MatrixXd R = MatrixXd(n_z, n_z);
    S = S + R_radar_;
    
    // calculate cross-correlation
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    
    for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd z_res = Zsig.col(i) - z_pred;
        NormalizeAngle(z_res, 1); // normalize phi if it is > |2*pi|
        VectorXd x_res = Xsig_pred_.col(i) - x_;
        NormalizeAngle(x_res, 3); // normalize phi if it is > |2*pi|
        Tc = Tc + weights_(i) * x_res * z_res.transpose();
    }
    
    // calculate Kalman gain
    MatrixXd K = Tc * S.inverse();
    
    // calculate innovation
    VectorXd innovation = z - z_pred;
    NormalizeAngle(innovation, 1); // normalize phi if it is > |2*pi|
    
    // update state mean and covariance
    x_ = x_ + K * innovation;
    P_ = P_ - K * S * K.transpose();
    
    // determine normalized innovation squared for testing parameters' goodness of fit
    NIS_radar_ = innovation.transpose() * S.inverse() * innovation;
}


/**
 * GenerateAugmentedSigmaPoints
 * @return Xsig MatrixXd with augmented sigma points
 */
MatrixXd UKF::GenerateAugmentedSigmaPoints(){
    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;
    
    // create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5, 5) = P_;
    P_aug(5, 5) = std_a_ * std_a_;
    P_aug(6, 6) = std_yawdd_ * std_yawdd_;
    
    // calculate square root of augmented state covariance
    MatrixXd L = P_aug.llt().matrixL();
    
    // create augmented sigma points
    MatrixXd Xsig(n_aug_, 2 * n_aug_ + 1);
    Xsig.col(0) = x_aug;
    
    for (unsigned int i = 0; i < n_aug_; i++)
    {
        Xsig.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }
    return Xsig;
}


/**
 * PredictAugmentedSigmaPoints
 * @param Xsig MatrixXd with augmented sigma points
 * @param dt Time between k and k+1 in s
 */
void UKF::PredictAugmentedSigmaPoints(MatrixXd Xsig, double dt){
    
    for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++){
        // write values from matrix to properly named variables for readability
        double px = Xsig(0, i);
        double py = Xsig(1, i);
        double v  = Xsig(2, i);
        double yaw = Xsig(3, i);
        double yawd = Xsig(4, i);
        double nu_a = Xsig(5, i);
        double nu_yawdd = Xsig(6, i);
        
        // predict state variables
        double px_p, py_p;
        
        if (fabs(yawd) > 0.001) {
            px_p = px + v / yawd * (sin(yaw + yawd * dt) - sin(yaw));
            py_p = py + v / yawd * (cos(yaw) - cos(yaw + yawd * dt));
        }
        else {
            px_p = px + v * dt * cos(yaw);
            py_p = py + v * dt * sin(yaw);
        }
        
        double v_p = v;
        double yaw_p = yaw + yawd * dt;
        double yawd_p = yawd;
        
        // add noise as described by Julier and Uhlmann (1997)
        // --> noise values evaluate to 0 (mean noise) or +/- sqrt(lambda_ + n_aug_) * noise standard deviation
        px_p = px_p + 0.5 * nu_a * dt * dt * cos(yaw);
        py_p = py_p + 0.5 * nu_a * dt * dt * sin(yaw);
        v_p = v_p + nu_a * dt;
        yaw_p = yaw_p + 0.5 * nu_yawdd * dt * dt;
        yawd_p = yawd_p + nu_yawdd * dt;
        
        // write values predicted state values to matrix
        Xsig_pred_(0, i) = px_p;
        Xsig_pred_(1, i) = py_p;
        Xsig_pred_(2, i) = v_p;
        Xsig_pred_(3, i) = yaw_p;
        Xsig_pred_(4, i) = yawd_p;
    }
}


/**
 * NormalizeAngle normalize angle if it is > |2*pi|
 * @param vector Vector containing relevant angle value at position idx
 * @param idx Position at which vector contains relevant angle value
 */
void UKF::NormalizeAngle(VectorXd vector, int idx){
    while (vector(idx) > M_PI) vector(idx) -= 2. * M_PI;
    while (vector(idx) < -M_PI) vector(idx) += 2. * M_PI;
}


