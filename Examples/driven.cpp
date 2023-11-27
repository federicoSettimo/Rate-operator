/*
  Ph cov in the x direction, driving sigma_z
*/
#include "../roqj.h"

using namespace std;
using namespace Eigen;

double gamma_p (double t) {return 1.;}
double gamma_m (double t) {return 1.;}
double b (double t) {return 10.;}

//Matrix2cd sigma_p_x{{-.5,.5},{-.5,.5}}, sigma_m_x = sigma_p_x.transpose();
Matrix2cd sigma_p_x = plus_state*minus_state.adjoint(), sigma_m_x = sigma_p_x.transpose();

MatrixXcd H (double t) {
  return b(t)*sigma_z;
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_p(t)*sigma_p_x*rho*sigma_m_x + gamma_m(t)*sigma_m_x*rho*sigma_p_x;
}

MatrixXcd Gamma (double t) {
  return gamma_p(t)*sigma_m_x*sigma_p_x + gamma_m(t)*sigma_p_x*sigma_m_x;
}

VectorXcd Phi (const VectorXcd &psi, double t, bool jumped) {
  complex<double> alpha = psi(1), phi_p, phi_m;
  double a = abs(alpha), sq = sqrt(1.-a*a), th = arg(alpha);

  if (jumped) {
    //cout << "Jumped\n";
    return 0.*ground_state;
  }


  //cout << "Not jumped, alpha = " << alpha << ", a = " << a << endl;
  //cout << "\tpsi = " << psi << endl;

  if (a >= .5) {
    phi_m = sq;
    phi_p = sq*(a*exp(2.*I*th) + a - exp(I*th)*conj(phi_m))/a;
  }
  else {
    phi_p = sq;
    phi_m = exp(-I*th)*a + a*exp(I*th)*(1. - conj(phi_p)/sq);
  }

  /*if (a >= .5) {
    phi_m = sq;
    phi_p = (sqrt(1. - pow(a,2))*(a + a*exp(2.*I*th) - exp(I*th)*conj(phi_m)))/a;
  }
  else {
    phi_p = sq;
    phi_m = (a*(1. + exp(2.*I*th) - (exp(2.*I*th)*conj(phi_p))/sqrt(1 - pow(a,2))))/exp(I*th);
  }*/

  return phi_m*ground_state + phi_p*excited_state;
}

double observable (const MatrixXcd &rho) {return real((rho*sigma_x).trace());}

int main () {
  double tmin = 0., tmax = 3., dt = 0.001;
  int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 10;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj);

  Vector2cd initialState = .5*excited_state + .1*ground_state;
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}