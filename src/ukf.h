/**
 * Unscented Kalman Filter as described in
 * Julier, S., and J. Uhlmann. 1997. A new extension of the Kalman filter to nonlinear systems
 */

#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"


class UKF {
    public:
        /**
        * Constructor
        */
        UKF();

    
        /**
        * Destructor
        */
        virtual ~UKF();

    
        /**
        * ProcessMeasurement
        * @param meas_package The latest measurement data of either radar or laser
        */
        void ProcessMeasurement(MeasurementPackage meas_package);

    
        /**
        * Prediction Predicts augmented sigma points, the state, and the state covariance
        * matrix
        * @param delta_t Time between k and k+1 in s
        */
        void Prediction(double delta_t);

    
        /**
        * Updates the state and the state covariance matrix using a laser measurement
        * @param meas_package The measurement at k+1
        */
        void UpdateLidar(MeasurementPackage meas_package);

    
        /**
        * Updates the state and the state covariance matrix using a radar measurement
        * @param meas_package The measurement at k+1
        */
        void UpdateRadar(MeasurementPackage meas_package);
    
    
        /**
         * GenerateAugmentedSigmaPoints
         * @param
         * @return MatrixXd with augmented sigma points
         */
        Eigen::MatrixXd GenerateAugmentedSigmaPoints();
    
    
        /**
         * PredictAugmentedSigmaPoints
         * @param Xsig MatrixXd with augmented sigma points
         * @param dt Time between k and k+1 in s
         */
        void PredictAugmentedSigmaPoints(Eigen::MatrixXd Xsig, double dt);

    
        /**
         * NormalizeAngle normalize angle if it is > |2*pi|
         * @param vector Vector containing relevant angle value at position idx
         * @param idx Position at which vector contains relevant angle value
         */
        void NormalizeAngle(Eigen::VectorXd vector, int idx);
    
        // initially set to false, set to true in first call of ProcessMeasurement
        bool is_initialized_;

        // if this is false, laser measurements will be ignored (except for init)
        bool use_laser_;

        // if this is false, radar measurements will be ignored (except for init)
        bool use_radar_;

        // state dimension
        int n_x_;

        // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
        Eigen::VectorXd x_;

        // state covariance matrix
        Eigen::MatrixXd P_;

        // augmented state dimension
        int n_aug_;

        // sigma point spreading parameter
        double lambda_;
    
        // weights of sigma points
        Eigen::VectorXd weights_;
    
        // predicted sigma points matrix
        Eigen::MatrixXd Xsig_pred_;

        // laser measurement noise standard deviation position1 in m
        double std_laspx_;

        // laser measurement noise standard deviation position2 in m
        double std_laspy_;

        // radar measurement noise standard deviation radius in m
        double std_radr_;

        // radar measurement noise standard deviation angle in rad
        double std_radphi_;

        // radar measurement noise standard deviation radius change in m/s
        double std_radrd_ ;
    
        // laser measurement noise covariance matrix
        Eigen::MatrixXd R_laser_;
    
        // radar measurement noise covariance matrix
        Eigen::MatrixXd R_radar_;

        // time when the state is true, in us
        long long previous_timestamp_;

        // process noise standard deviation longitudinal acceleration in m/s^2
        double std_a_;

        // process noise standard deviation yaw acceleration in rad/s^2
        double std_yawdd_;
    
        // current normalized innovation squared (NIS) of laser
        double NIS_laser_;
    
        // current normalized innovation squared (NIS) of radar
        double NIS_radar_;
};

#endif  // UKF_H
