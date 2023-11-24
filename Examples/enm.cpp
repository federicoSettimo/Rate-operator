// Eternally non-Markovian dynamcis
#include "../roqj.h"

using namespace std;
using namespace Eigen;

double gamma_p (double t) {return 1.;}
double gamma_m (double t) {return 1.;}
double gamma_z (double t) {return -.5*tanh(t);}

MatrixXcd H (double t) {
  return 0.*sigma_x;
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_p(t)*sigma_p*rho*sigma_m + gamma_m(t)*sigma_m*rho*sigma_p + gamma_z(t)*sigma_z*rho*sigma_z;
}

MatrixXcd Gamma (double t) {
  return gamma_p(t)*sigma_m*sigma_p + gamma_m(t)*sigma_p*sigma_m + gamma_z(t)*id;
}

VectorXcd Phi (const VectorXcd &psi, double t, bool jumped) {
  double a = real(psi(1)), gm = gamma_m(t), gp = gamma_p(t), gz = gamma_z(t);

  double phi1ub = pow((1.-a*a),1.5)/(a*a)*gm + 3.*sqrt(1.-a*a)*gz;
  double phi1lb = -a*a/sqrt(1.-a*a)*gp - gz*sqrt(1.-a*a);

  // Use either one of the following two lines to have jumps either only to |0> or only to |1>
  // Works fine with any phi = lambda*phi1ub + (1.-lambda)*phi1lb, for 0 < lambda < 1
  //double phi1 = phi1ub;
  double phi1 = phi1lb;

  if (!jumped)
    return a*(2.*gz - phi1/sqrt(1.-a*a))*ground_state + phi1*excited_state;
  return -gz*psi;
}

double observable (const MatrixXcd &rho) {return real((rho*sigma_z).trace());}

int main () {
  double tmin = 0., tmax = 2, dt = 0.01;
  int N_ensemble = 1000, Ncopies = 3, dimH = 2, Ntraj = 7;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj, true, 1e-2);

  Vector2cd initialState = .9*excited_state - .5*ground_state;
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}