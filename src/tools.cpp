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
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd res(4);
  res << 0, 0, 0, 0;

  if(estimations.size() == 0 || estimations.size() != ground_truth.size()) {
    cerr << "Wrong estimations size or length...\n";
    return res;
  }

  for(int i=0; i<estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()  * residual.array();
    res += residual;
  }
  res /= (double)estimations.size();
  res = res.array().sqrt();

  return res;
}

MatrixXd Tools::DiagonizeMatrix(const MatrixXd& mat)
{
  MatrixXd m = MatrixXd::Zero(mat.rows(), mat.cols());
  m.diagonal() = mat.diagonal();

  return m;
}

double Tools::NormalizeAngle(const double angle)
{
  //std::cout << "Angle Normalization:\n\tbefore: " << angle;
  //double a = angle;
  //while (a > M_PI) a -= 2.*M_PI;
  //while (a < -M_PI) a += 2.*M_PI;
  double a = atan2(sin(angle), cos(angle));
  //std::cout << "\n\tafter: " << a << std::endl;
  return a;
}
