#ifndef DOUBLEMASS_H
#define DOUBLEMASS_H

#include "Equation.h"

class DoubleMass : public Equation {
public:
    DoubleMass();
    ~DoubleMass();

    virtual void RHS(double t, const double *q, double *dq) const;

private:
    double &A_, &B_, &J_, &mp_, &r_, &R_, &Omega_, &Gamma0_, &epsGamma_, &delta_, &mu1_, &mu2_, &mu3_;
    mutable double K_[3][3], invK_[3][3];
};

#endif
