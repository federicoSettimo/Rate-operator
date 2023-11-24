// Phase covariant violating P divisibility from gamma_z < sqrt(gamma_p*gamma_m)/2
#include "../roqj.h"

using namespace std;
using namespace Eigen;

double gamma_p (double t) {return exp(-.5*t);}
double gamma_m (double t) {return exp(-.25*t);}
double gamma_z (double t) {return 4*exp(-3.*t/8.)*cos(2.*t)*.5;}

//double b (double t) {return 0.5*(1. + erf(2.*sqrt(2.)*(t-1.)));}
double b (double t) {return 0.;}

MatrixXcd H (double t) {
  return -0.5*b(t)*sigma_z;
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_p(t)*sigma_p*rho*sigma_m + gamma_m(t)*sigma_m*rho*sigma_p + gamma_z(t)*sigma_z*rho*sigma_z;
}

MatrixXcd Gamma (double t) {
  return gamma_p(t)*sigma_m*sigma_p + gamma_m(t)*sigma_p*sigma_m + gamma_z(t)*id;
}

MatrixXcd C (const VectorXcd &psi, double t) {
  double mu, c3, a2 = norm(psi(1)), eps = -.5*sqrt(gamma_m(t)*gamma_p(t)) - gamma_z(t);
  eps = eps > 0. ? eps : 0.;
  mu = a2 == 0. ? sqrt(gamma_p(t)*gamma_m(t)) + 2.*eps : 2.*gamma_z(t) + gamma_m(t)*(1-a2)/a2;

  mu = a2 < .5 ? ((1.-a2)*sqrt(gamma_p(t)*gamma_m(t)) - a2*gamma_p(t))/(1.-a2) + 2.*eps : 2.*gamma_z(t) + gamma_m(t)*(1-a2)/a2;

  c3 = gamma_z(t) - mu;
  return 2.*mu*sigma_p*sigma_m + c3*id;
}

VectorXcd Phi (const VectorXcd &psi, double t, bool jumped) {return C(psi,t)*psi;}

double observable (const MatrixXcd &rho) {return real((rho*sigma_z).trace());}

int main () {
  double tmin = 0., tmax = 3., dt = 0.01;
  int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 5;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj);
  Vector2cd initialState = cos(M_PI/8.)*ground_state + sin(M_PI/8.)*excited_state;
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}