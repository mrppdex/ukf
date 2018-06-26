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
  std_a_ = 0.4;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 50;

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

  lambda_ = 3 - n_aug_;

  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
  0, std_radphi_*std_radphi_, 0,
  0, 0, std_radrd_*std_radrd_;

  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;



  //set vector for weights
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_[0] = weight_0;
  for (int i=1; i< 2*n_aug_ + 1; ++i) {
      double weight = 0.5/(n_aug_ + lambda_);
      weights_[i] = weight;
  }

/*
  P_ = MatrixXd::Identity(n_x_, n_x_).array()*1000;
  P_(0, 0) = 1;
  P_(1, 1) = 1;
  P_(n_x_ - 1, n_x_ - 1) = 10;
*/

  P_ = MatrixXd(n_x_, n_x_);
  P_ <<   0.03,  0,  0,  0,  0,
          0,  0.03,  0,  0,  0,
          0,  0, 1,  0,  0,
          0,  0,  0, 1 ,  0,
          0,  0,  0,  0,20;

  Q_ = MatrixXd::Zero(2, 2);
  Q_(0, 0) = std_a_*std_a_;
  Q_(1, 1) = std_yawdd_*std_yawdd_;

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
  std::cout << "Measured x:\n" << meas_package.raw_measurements_ << std::endl;
  if (!is_initialized_) {
    x_ = VectorXd(5);
    x_ << 0, 0, 5, 0, -180;
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      x_(0) = rho*cos(phi);
      x_(1) = rho*sin(phi);
      //x_(2) = 3;
    }
    is_initialized_ = true;
    previous_timestamp_ = meas_package.timestamp_;
    return;
  }

  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1.e6;
  Prediction(dt);

  if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    std::cout << "UpdateLidar\n";
    UpdateLidar(meas_package);
  } else if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    std::cout << "UpdateRadar\n";
    UpdateRadar(meas_package);
  }

  //P_ = tools.DiagonizeMatrix(P_);


  std::cout << "-----------------------------------------------------------\n";
  std::cout << P_ << std::endl;
  std::cout << "-----------------------------------------------------------\n";

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
  std::cout << "Predicted x:\n" << x_ << std::endl;
  std::cout << "GenerateSigmaPoints\n";
  MatrixXd X_sigma = GenerateSigmaPoints();
  std::cout << "PredictSigmaPoints\n";
  Xsig_pred_ = PredictSigmaPoints(X_sigma, delta_t);
  std::cout << "PredictMeanAndCovariance\n";
  PredictMeanAndCovariance(Xsig_pred_, &x_, &P_);

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  VectorXd z = meas_package.raw_measurements_;
  VectorXd y = z - H_*x_;
  std::cout << "\t\ty\n";
  MatrixXd S = H_*P_*H_.transpose() + R_lidar_; //2x5*5x5*5x2 + 2x2;
  std::cout << "\t\tS\n";
  MatrixXd K = P_*H_.transpose()*S.inverse();
  std::cout << "\t\tK\n";
  x_ = x_ + K*y;
  std::cout << "\t\tx_ again\n";
  P_ = (MatrixXd::Identity(n_x_, n_x_) - K*H_)*P_;
  std::cout << "\t\tP_\n";
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
  VectorXd z_measured = meas_package.raw_measurements_;
  std::cout << "\tPredictRadarMeasurement\n";
  Zsig_ = PredictRadarMeasurement(Xsig_pred_, &z_, &S_);
  std::cout << "\tUpdateRadarState\n";
  UpdateRadarState(Xsig_pred_, z_measured, &x_, &P_);
}

MatrixXd UKF::GenerateSigmaPoints()
{
  MatrixXd X_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  MatrixXd P_aug = MatrixXd::Zero(n_x_ + 2, n_x_ + 2);

  P_aug << P_, MatrixXd::Zero(n_x_, 2),
           MatrixXd::Zero(2, n_x_), Q_;
  std::cout << "\t P_aug initialized...\n";
  MatrixXd P_l = P_aug.llt().matrixL();

  VectorXd x_aug_ = VectorXd::Zero(n_aug_);
  x_aug_ << x_, 0, 0;

  X_aug.col(0) = x_aug_; //k|k

  for(int i=0; i < n_aug_; ++i) {
    X_aug.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_)*P_l.col(i);
    X_aug.col(i + n_aug_ + 1) = x_aug_ - sqrt(lambda_ + n_aug_)*P_l.col(i);
  }

  return X_aug;
}

VectorXd UKF::PredictState(const VectorXd &x, const double delta_t)
{
  const double px = x[0];
  const double py = x[1];
  const double v  = x[2];
  const double phi = x[3];
  const double phi_dot = x[4];
  const double nu_a = x[5];
  const double nu_phi = x[6];

  VectorXd x_tmp = VectorXd(n_x_);

  VectorXd nu_k = VectorXd(n_x_);
  nu_k << 0.5*delta_t*delta_t*cos(phi)*nu_a,
  0.5*delta_t*delta_t*sin(phi)*nu_a,
  delta_t*nu_a,
  0.5*delta_t*delta_t*nu_phi,
  delta_t*nu_phi;

  if(phi_dot > 0.0001) {
      x_tmp[0] = px + v/phi_dot*(sin(phi+phi_dot*delta_t) - sin(phi));
      x_tmp[1] = py + v/phi_dot*(-cos(phi+phi_dot*delta_t) + cos(phi));
      x_tmp[2] = v + 0;
      x_tmp[3] = phi + phi_dot*delta_t;
      x_tmp[4] = phi_dot + 0;
  } else {
      x_tmp[0] = px + v*cos(phi)*delta_t;
      x_tmp[1] = py + v*sin(phi)*delta_t;
      x_tmp[2] = v + 0;
      x_tmp[3] = phi + 0;
      x_tmp[4] = phi_dot + 0;
  }

  x_tmp += nu_k;
  return x_tmp;
}

MatrixXd UKF::PredictSigmaPoints(const MatrixXd &Xsig_aug, const double delta_t)
{
    MatrixXd Xsig_pred = MatrixXd(n_x_, 2*n_aug_ + 1);

    for(int i=0; i<Xsig_aug.cols(); ++i) {
        VectorXd x = Xsig_aug.col(i);
        Xsig_pred.col(i) = PredictState(x, delta_t);
    }

    return Xsig_pred;
}


void UKF::PredictMeanAndCovariance(const MatrixXd &Xsig_pred, VectorXd* x_out, MatrixXd* P_out)
{

    //create vector for weights
    //VectorXd weights = VectorXd(2*n_aug_+1);

    //create vector for predicted state
    VectorXd x = VectorXd::Zero(n_x_); //VectorXd(n_x);

    //create covariance matrix for prediction
    MatrixXd P = MatrixXd::Zero(n_x_, n_x_); //MatrixXd(n_x, n_x);


    //predict state mean
    for(int i=0; i<Xsig_pred.cols(); ++i) {
        x += weights_[i]*Xsig_pred.col(i);
    }
    std::cout << "\tstate mean calculated...\n";

    //predict state covariance matrix
    for(int i=0; i<Xsig_pred.cols(); ++i) {
        VectorXd v_mean = Xsig_pred.col(i) - x;

        //normalize yaw
        v_mean(3) = atan2(sin(v_mean(3)), cos(v_mean(3)));

        MatrixXd v_p = v_mean*v_mean.transpose();
        v_p *= weights_[i];
        P += v_p;

    }
    std::cout << "\tcovariance matrix calculated...\n";

    *x_out = x;
    *P_out = P;
    return;
}

MatrixXd UKF::PredictRadarMeasurement(const MatrixXd &Xsig_pred, VectorXd* z_out, MatrixXd* S_out)
{
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

    //mean predicted measurement
    VectorXd z_pred = VectorXd::Zero(n_z_);

    //measurement covariance matrix S
    MatrixXd S = MatrixXd::Zero(n_z_,n_z_);

    //transform sigma points into measurement space
    //calculate mean predicted measurement
    //calculate innovation covariance matrix S
    for(int i=0; i<Xsig_pred.cols(); ++i) {
        const double px = Xsig_pred(0, i);
        const double py = Xsig_pred(1, i);
        const double v = Xsig_pred(2, i);
        const double yaw = Xsig_pred(3, i);
        //const double yawd = Xsig_pred(4, i);

        VectorXd z_tmp = VectorXd(z_pred.size());
        z_tmp[0] = sqrt(px*px + py*py);
        z_tmp[1] = atan2(py, px);
        if(fabs(z_tmp[0]) > 0.001) {
            z_tmp[2] = (px*v*cos(yaw) + py*v*sin(yaw))/z_tmp[0];
        } else {
            z_tmp[2] = 0.0;
        }

        Zsig.col(i) = z_tmp;
        z_pred += weights_[i]*z_tmp;
    }

    for(int i=0; i<Zsig.cols(); ++i) {
        VectorXd v = Zsig.col(i);
        MatrixXd v_mean = (v - z_pred)*(v - z_pred).transpose();
        v_mean *= weights_[i];
        S += v_mean;
    }

    S += R_radar_;

    *z_out = z_pred;
    *S_out = S;

    return Zsig;
}

void UKF::UpdateRadarState(const MatrixXd &Xsig_pred, const VectorXd z, VectorXd* x_out, MatrixXd* P_out)
{

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_);

    //calculate cross correlation matrix
    //calculate Kalman gain K;
    //update state mean and covariance matrix
    for(int i=0; i<Zsig_.cols(); ++i) {
        Tc += weights_(i)*(Xsig_pred_.col(i) - x_)*(Zsig_.col(i) - z_).transpose();
    }
    std::cout << "\t\tfirst loop finished...\n";

    MatrixXd K = Tc*S_.inverse();

    x_ = x_ + K*(z - z_);
    x_(3) = atan2(sin(x_(3)), cos(x_(3)));

    std::cout << "\t\tx_ updated\n";
    P_ = P_ - K*S_*K.transpose();
    std::cout << "\t\tP_ updated\n";


}
