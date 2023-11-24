// Non-Markovianity from gamma_p = gamma_m < 0
// Writing everyting in the |+-> basis
#include "../roqj.h"

using namespace std;
using namespace Eigen;

Matrix2cd sigma_p_pm {{1,-1},{1,-1}}, sigma_m_pm {{1,1},{-1,-1}}, sigma_z_pm {{0,1},{1,0}}, sigma_x_pm {{1,0},{0,-1}}, sigma_y_pm {{0,I},{-I,0}};
Vector2cd plus_pm = excited_state, minus_pm = ground_state, excited_state_pm = (plus_pm+minus_pm)/sqrt(2.), ground_state_pm = (plus_pm-minus_pm)/sqrt(2.);

// Note: it must be gp + gm + gz >= 0 at all times
double kappa = .25;
double gamma_p (double t) {return .5*exp(-.1*t)*(kappa + (1.-kappa)*exp(-.25*t)*cos(2.*t));}
double gamma_m (double t) {return gamma_p(t);}
double gamma_z (double t) {return .5*1.;}

MatrixXcd H (double t) {
  return MatrixXcd::Zero(2,2);
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_p(t)*sigma_p_pm*rho*sigma_m_pm + gamma_m(t)*sigma_m_pm*rho*sigma_p_pm + gamma_z(t)*sigma_z_pm*rho*sigma_z_pm;
}

MatrixXcd Gamma (double t) {
  return gamma_p(t)*sigma_m_pm*sigma_p_pm + gamma_m(t)*sigma_p_pm*sigma_m_pm + gamma_z(t)*id;
}

MatrixXcd C (const VectorXcd &psi, double t) {
  Vector2cd cpsi = psi;
  if (real(psi(0)) < 0.) cpsi -= psi;
  double a = real(cpsi(1)), cp, cm; // a = <psi|->
  Matrix2cd CC = MatrixXcd::Zero(2,2);
  double a2 = norm(a), g = gamma_m(t), gz = gamma_z(t);

  double cub, clb, c;

  clb = (2.*g + a2*gz)/(-1 + a2);
  cub = (2.*g + 8.*a2*g + gz - 3.*a2*gz)/a2;
  if (a2 <= .01) c = clb;
  else if (a2 >= .99) c = cub;
  else  if (a > 1./sqrt(2.)) c = cub;
  else c = clb;

  if (a2 >= .99) c = cub;
  else c = clb;

  CC(0,0) = c;
  CC(1,1) = 8.*g - 2.*gz - c;

  return CC;
}

VectorXcd Phi (const VectorXcd &psi, double t, bool jumped) {return C(psi,t)*psi;}

double observable (const MatrixXcd &rho) {return real((rho*sigma_z).trace());}

int main () {
  double tmin = 0., tmax = 3., dt = 0.01;
  int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 5;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj);
  Vector2cd initialState = cos(M_PI/8.)*ground_state_pm + sin(M_PI/8.)*excited_state_pm;
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}