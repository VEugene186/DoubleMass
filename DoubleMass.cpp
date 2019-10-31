#include "DoubleMass.h"
#include <cmath>

DoubleMass::DoubleMass() : Equation(6, 13),
        A_(parameters_[0]), B_(parameters_[1]), J_(parameters_[2]), mp_(parameters_[3]),
        r_(parameters_[4]), R_(parameters_[5]), Omega_(parameters_[6]), Gamma0_(parameters_[7]),
        epsGamma_(parameters_[8]), delta_(parameters_[9]), mu1_(parameters_[10]), mu2_(parameters_[11]),
        mu3_(parameters_[12]) 
{

}

DoubleMass::~DoubleMass() {

}

void DoubleMass::RHS(double t, const double *q, double *dq) const {
    double p1 = q[0];
    double p2 = q[1];
    double M = q[2];
    double phi = q[5];

    double c = cos(2.0 * M_PI * t);

    K_[0][0] = A_;       K_[0][1] = 0.0; K_[0][2] = 2 * mp_ * r_ * c;

    K_[1][0] = 0.0;      K_[1][1] = B_;  K_[1][2] = 0.0;

    K_[2][0] = K_[0][2]; K_[2][1] = 0.0; K_[2][2] = J_;

    double det = K_[1][1] * (K_[0][0] * K_[2][2] - K_[0][2] * K_[2][0]);
    invK_[0][0] = K_[1][1] * K_[2][2];
    invK_[0][1] = 0.0;
    invK_[0][2] = - K_[1][1] * K_[2][0];

    invK_[1][0] = 0.0;
    invK_[1][1] = K_[0][0] * K_[2][2] - K_[0][2] * K_[2][0];
    invK_[1][2] = 0.0;

    invK_[2][0] = invK_[0][2];
    invK_[2][1] = 0.0;
    invK_[2][2] = K_[0][0] * K_[1][1];
    
    double v1 = (invK_[0][0] * p1 + invK_[0][2] * (M + 2.0 * mp_ * r_ * R_ * Omega_ * c)) / det;
    double v2 = invK_[1][1] * p2 / det;
    double omega = (invK_[2][0] * p1 + invK_[2][2] * (M + 2.0 * mp_ * r_ * R_ * Omega_ * c)) / det;

    double Gamma = Gamma0_ + epsGamma_ * sin(2.0 * M_PI * t + delta_);

    double k = 2.0 * M_PI / Omega_;

    double sP, cP;
    sincos(phi, &sP, &cP);
    dq[0] = k * (p2 * omega - Gamma * v2 - mu1_ * v1);
    dq[1] = k * (-p1 * omega + Gamma * v1 - mu2_ * v2);
    dq[2] = k * (p1 * v2 - p2 * v1 - mu3_ * omega);
    dq[3] = k * (v1 * cP - v2 * sP);
    dq[4] = k * (v1 * sP + v2 * cP);
    dq[5] = k * omega;


}
