#include "BifurcationTree.h"
#include <cstdio>

BifurcationTree::BifurcationTree() {

}

BifurcationTree::~BifurcationTree() {

}

void BifurcationTree::calculate(Equation * eqs, int parNum, double begin, double end, int N, const double * initialPoint, int preIters, int saveCount) {
    int dim = eqs->getDim();
    method_.init(eqs);
    parValues_.resize(N);
    points_.resize(N);

    double step = (end - begin) / (N - 1);
    double * q0 = new double[dim];
    double * q1 = new double[dim];
    for (int j = 0; j < dim; j++) {
        q0[j] = initialPoint[j];
    }

    for (int i = 0; i < N; i++) {
        double value = begin + i * step;
        printf("it = %6d | val = %lg\n", i, value);

        eqs->setParameter(parNum, value);
        for (int k = 0; k < preIters; k++) {
            method_.map(eqs, 0.0, q0, q1, 10, 1.0);
            swap(q0, q1);
        }
        
        add(value, q0, dim);
        for (int k = 0; k < saveCount; k++) {
            method_.map(eqs, 0.0, q0, q1, 10, 1.0);
            swap(q0, q1);
            add(value, q0, dim);
        }
    }

    delete [] q0;
    delete [] q1;
}

void BifurcationTree::add(double value, const double * q, int dim) {
    parValues_.push_back(value);
    valarray<double> tmp(dim);
    for (int i = 0; i < dim; i++) {
        tmp[i] = q[i];
    }
    points_.push_back(tmp);
}

void BifurcationTree::saveToFile(const char * fileName) {
    int N = (int)parValues_.size();
    int dim = (int)points_[0].size();
    FILE * f = fopen(fileName, "w");
    for (int n = 0; n < N; n++) {
        fprintf(f, "%.15lg", parValues_[n]);
        for (int i = 0; i < dim; i++) {
            fprintf(f, "\t%.15lg", points_[n][i]);
        }
        fprintf(f, "\n"); 
    }
    fclose(f);
}
