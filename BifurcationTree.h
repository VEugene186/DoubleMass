#ifndef BIFURCATIONTREE_H
#define BIFURCATIONTREE_H

#include "Equation.h"
#include "RungeKutta.h"
#include <valarray>
#include <vector>

using namespace std;

class BifurcationTree {
public:
    BifurcationTree();
    ~BifurcationTree();

    void calculate(Equation * eqs, int parNum, double begin, double end, int N, const double * initialPoint, int preIters, int saveCount);
    void saveToFile(const char * fileName);
private:
    void add(double value, const double * q, int dim);
private:
    RungeKutta method_;
    vector<double> parValues_;
    vector<valarray<double > > points_;
};

#endif
