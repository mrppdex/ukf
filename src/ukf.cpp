#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  tools = Tools();

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:
  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */

  n_z_ = 3;
  n_x_ = 5;
  n_aug_ = 7;

  lambda_ = 3 - n_aug_; //3

  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_,     0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,   std_radrd_*std_radrd_;

  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;



  //set vector for weights
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i< 2*n_aug_ + 1; ++i) {
      double weight = 0.5/(n_aug_ + lambda_);
      weights_(i) = weight;
  }

  //object covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ <<   0.03,  0,  0,  0,  0,
          0,  0.03,  0,  0,  0,
          0,     0,  1,  0,  0,
          0,     0,  0,  1 , 0,
          0,     0,  0,  0 , 1;

  //process covariance matrix
  Q_ = MatrixXd::Zero(2, 2);
  Q_(0, 0) = std_a_*std_a_;
  Q_(1, 1) = std_yawdd_*std_yawdd_;

  //lidar measurement matrix
  H_ = MatrixXd(2, n_x_);
  H_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;

}

UKF::~UKF() {}

/**
 * @param meas_package meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_) {
    //initialize state vector

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      x_(0) = rho*cos(phi);
      x_(1) = rho*sin(phi);
    }
    is_initialized_ = true;
    previous_timestamp_ = meas_package.timestamp_;
    return;
  }

  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1.e6;
  Prediction(dt);

  if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  } else if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }

  NormalizedInnovationSquared(meas_package.raw_measurements_);

  //diagonalize object matrix
  //P_ = tools.DiagonizeMatrix(P_);

  previous_timestamp_ = meas_package.timestamp_;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */


  MatrixXd X_sigma_aug = GenerateSigmaPoints();
  Xsig_pred_ = PredictSigmaPoints(X_sigma_aug, delta_t);
  PredictMeanAndCovariance();

}

MatrixXd UKF::GenerateSigmaPoints()
{
  MatrixXd X_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  MatrixXd P_aug = MatrixXd::Zero(n_x_ + 2, n_x_ + 2);

  //create augmented covariance matrix
  P_aug << P_, MatrixXd::Zero(n_x_, 2),
           MatrixXd::Zero(2, n_x_), Q_;

  //square root of P (cholesky)
  MatrixXd P_l = P_aug.llt().matrixL();

  VectorXd x_aug_ = VectorXd::Zero(n_aug_);
  x_aug_ << x_, 0, 0;

  X_aug.col(0) = x_aug_; //k|k

  // find sigma points
  for(int i=0; i < n_aug_; ++i) {
    X_aug.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_)*P_l.col(i);
    X_aug.col(i + n_aug_ + 1) = x_aug_ - sqrt(lambda_ + n_aug_)*P_l.col(i);
  }

  return X_aug;
}

MatrixXd UKF::PredictSigmaPoints(const MatrixXd &Xsig_aug, const double delta_t)
{
    /**
     * Predicts sigma points for augmented State
     */

    MatrixXd Xsig_pred = MatrixXd(n_x_, 2*n_aug_ + 1);

    for(int i=0; i<Xsig_aug.cols(); ++i) {
        Xsig_pred.col(i) = PredictState(Xsig_aug.col(i), delta_t);
    }

    return Xsig_pred;
}

VectorXd UKF::PredictState(const VectorXd &x_aug_vec, const double delta_t)
{
  const double px = x_aug_vec[0];
  const double py = x_aug_vec[1];
  const double v  = x_aug_vec[2];
  const double yaw = x_aug_vec[3];
  const double yawd = x_aug_vec[4];
  const double nu_a = x_aug_vec[5];
  const double nu_yawdd = x_aug_vec[6];

  VectorXd x_tmp = VectorXd(n_x_);

  VectorXd nu_k = VectorXd(n_x_);

  //nu_k {px, py, v, yaw, yaw_rate}
  nu_k << 0.5*delta_t*delta_t*cos(yaw)*nu_a,
          0.5*delta_t*delta_t*sin(yaw)*nu_a,
          delta_t*nu_a,
          0.5*delta_t*delta_t*nu_yawdd,
          delta_t*nu_yawdd;

  //update x, avoid division by 0
  if(fabs(yawd) > 1.e-6) {
      x_tmp[0] = px + v/yawd*(sin(yaw+yawd*delta_t) - sin(yaw));
      x_tmp[1] = py + v/yawd*(-cos(yaw+yawd*delta_t) + cos(yaw));
      x_tmp[2] = v + 0;
      x_tmp[3] = yaw + yawd*delta_t;
      x_tmp[4] = yawd + 0;
  } else {
      x_tmp[0] = px + v*cos(yaw)*delta_t;
      x_tmp[1] = py + v*sin(yaw)*delta_t;
      x_tmp[2] = v + 0;
      x_tmp[3] = yaw + 0;
      x_tmp[4] = yawd + 0;
  }

  x_tmp += nu_k;
  return x_tmp;
}


void UKF::PredictMeanAndCovariance()
{

    //create vector for weights
    //VectorXd weights = VectorXd(2*n_aug_+1);

    //create vector for predicted state
    VectorXd x = VectorXd::Zero(n_x_); //VectorXd(n_x);

    //create covariance matrix for prediction
    MatrixXd P = MatrixXd::Zero(n_x_, n_x_); //MatrixXd(n_x, n_x);


    //predict state mean
    for(int i=0; i<Xsig_pred_.cols(); ++i) {
        x += weights_(i)*Xsig_pred_.col(i);
    }


    //predict state covariance matrix
    for(int i=0; i<Xsig_pred_.cols(); ++i) {
        VectorXd x_diff = Xsig_pred_.col(i) - x;

        //normalize yaw
        x_diff(3) = tools.NormalizeAngle(x_diff(3));
        P += weights_(i)*x_diff*x_diff.transpose();

    }

    x_ = x;
    P_ = P;
    return;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
  Kalman filter algorithm
  */

  VectorXd z = meas_package.raw_measurements_;
  VectorXd y = z - H_*x_;

  S_ = H_*P_*H_.transpose() + R_lidar_; //dimensions 2x5*5x5*5x2 + 2x2 = 2x2;
  MatrixXd K = P_*H_.transpose()*S_.inverse();
  x_ = x_ + K*y;
  P_ = (MatrixXd::Identity(n_x_, n_x_) - K*H_)*P_;

  return;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /**
  TODO:
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the radar NIS.
  */
  Zsig_ = PredictRadarMeasurement();
  UpdateRadarState(meas_package.raw_measurements_);
}

MatrixXd UKF::PredictRadarMeasurement()
{
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

    //transform sigma points into measurement space
    //calculate mean predicted measurement
    //calculate innovation covariance matrix S
    for(int i=0; i < 2*n_aug_ + 1; ++i) {
        const double px = Xsig_pred_(0, i);
        const double py = Xsig_pred_(1, i);
        const double v = Xsig_pred_(2, i);
        const double yaw = Xsig_pred_(3, i);
        //const double yawd = Xsig_pred(4, i);

        VectorXd z_tmp = VectorXd(n_z_);
        z_tmp[0] = sqrt(px*px + py*py);
        z_tmp[1] = atan2(py, px);
        //avoid division by zero
        if(fabs(z_tmp[0]) > 0.001) {
            z_tmp[2] = (px*v*cos(yaw) + py*v*sin(yaw))/z_tmp[0];
        } else {
            z_tmp[2] = 0.0;
        }

        Zsig.col(i) = z_tmp;
        //z_pred += weights_[i]*z_tmp;
    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd::Zero(n_z_);
    for (int i=0; i < 2*n_aug_ + 1; ++i) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd::Zero(n_z_,n_z_);
    for(int i=0; i< 2*n_aug_ + 1; ++i) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        z_diff(1) = tools.NormalizeAngle(z_diff(1));
        S += weights_(i) * z_diff * z_diff.transpose();
    }

    S += R_radar_;

    z_  = z_pred;
    S_ = S;

    return Zsig;
}

void UKF::UpdateRadarState(const VectorXd& z)
{
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_);

  //calculate cross correlation matrix
  //calculate Kalman gain K;
  //update state mean and covariance matrix

  for(int i=0; i < 2*n_aug_ + 1; ++i) {
    VectorXd z_diff = Zsig_.col(i) - z_;
    z_diff(1) = tools.NormalizeAngle(z_diff(1));

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = tools.NormalizeAngle(x_diff(3));

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }


  MatrixXd K = Tc * S_.inverse();
  VectorXd z_diff = z - z_;
  z_diff(1) = tools.NormalizeAngle(z_diff(1));

  x_ = x_ + K*z_diff;

  P_ = P_ - K*S_*K.transpose();

}

void UKF::NormalizedInnovationSquared(const VectorXd& x)
{
  //chi-squared threshold for 95% = 7.8
  const double kThreshold = {7.8};
  unsigned long int counter = {0};
  double epsilon = {0.0};

  //check if got lidar or radar measurement and calculate epsilon
  if (use_laser_ && x.size() == 2) {
    epsilon = (x - x_.head(2)).transpose() * S_.inverse() * (x - x_.head(2));
  } else if (use_radar_ && x.size() == 3) {
    epsilon = (x - z_).transpose() * S_.inverse() * (x - z_);
  }
  //std::cout << "epsilon = " << epsilon << std::endl;
  nis.push_back(epsilon);

  for(auto a: nis) {
    if (a < kThreshold) ++counter;
  }

  std::cout << "NIS < 7.8 = " << counter/(double)nis.size() << std::endl;

}
