/*
  Ph cov in the x direction, driving sigma_z
*/
#include "../roqj.h"

using namespace std;
using namespace Eigen;

double sgn (double x) {return x >= 0. ? 1. : 0.;}

double gamma_p (double t) {return 1.;}
double gamma_m (double t) {return 1.;}
double gamma_x (double t) {return -.5*tanh(t);}
double b (double t) {return 10.;}

Matrix2cd sigma_p_x = plus_state*minus_state.adjoint(), sigma_m_x = sigma_p_x.transpose();

MatrixXcd H (double t) {
  return b(t)*sigma_z;
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_p(t)*sigma_p_x*rho*sigma_m_x + gamma_m(t)*sigma_m_x*rho*sigma_p_x + gamma_x(t)*sigma_x*rho*sigma_x;
}

MatrixXcd Gamma (double t) {
  return gamma_p(t)*sigma_m_x*sigma_p_x + gamma_m(t)*sigma_p_x*sigma_m_x + gamma_x(t)*sigma_x*sigma_x;
}

VectorXcd Phi (const VectorXcd &psi, double t, bool jumped) {
  complex<double> alpha = psi(1), phi_e, phi_g;
  double a = abs(alpha), sq = sqrt(1.-a*a), th = arg(alpha), gx = gamma_x(t);

  if (jumped) {
    return 0.*ground_state;
  }

  double xi_lb = -(1.+2.*gx)*(1.-a*a) - a*a;
  double xi_ub = real((1.-2.*gx)*(conj(alpha)*conj(alpha) + alpha*alpha) + 3.*a*a + (a*a*a*a)/(1.-a*a)*(1.+2.*gx));

  //double xi = xi_lb; // Gives dynamics with |psi_det> spiraling towards |1>
  //double xi = xi_ub; // Gives dynamics with |psi_det> spiraling towards |0>
  double xi = .5*(xi_lb+xi_ub); // |psi_det> remains at the same z component as rho asymptotically

  phi_g = xi*exp(I*th)/(2.*a);
  phi_e = (sqrt(1 - pow(abs(alpha),2))*(alpha - 2.*alpha*gx + conj(alpha) - conj(phi_g)))/conj(alpha);

  return phi_g*ground_state + phi_e*excited_state;
}

double observable (const MatrixXcd &rho) {return real((rho*sigma_x).trace());}

int main () {
  double tmin = 0., tmax = 3., dt = 0.001;
  int N_ensemble = 1000, Ncopies = 3, dimH = 2, Ntraj = 10;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj);

  Vector2cd initialState = .5*excited_state + .1*ground_state;
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}