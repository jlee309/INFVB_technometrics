#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::MatrixXd myfun(Eigen::Map<Eigen::MatrixXd> a,
                      Eigen::Map<Eigen::MatrixXd> b) {
  Eigen::MatrixXd res = a * b;
  return res;
}

// [[Rcpp::export]]
Eigen::MatrixXd myfun2(Eigen::Map<Eigen::MatrixXd> a,
                      Eigen::Map<Eigen::MatrixXd> b) {
  Eigen::MatrixXd res = a.adjoint() * b;
  return res;
}

// [[Rcpp::export]]
Eigen::VectorXd myfun3(Eigen::Map<Eigen::MatrixXd> a,
                       Eigen::Map<Eigen::MatrixXd> b) {
  Eigen::MatrixXd int1 = a.adjoint() * b;
  Eigen::MatrixXd int2 = int1.array() * b.array();
  Eigen::VectorXd  res=int2.colwise().sum();
  return res;
}
