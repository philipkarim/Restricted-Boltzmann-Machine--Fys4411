#pragma once
#include "system.h"
#include <armadillo>

using namespace arma;

class SGD {
public:
    SGD(System* system, double eta);
    void SGDOptimize(arma::vec parameters_derivative);

private:
    double m_eta;
    vec a;
    vec b;
    mat w;

    class System* m_system = nullptr;
};