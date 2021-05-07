#pragma once
#include "system.h"
#include <armadillo>

using namespace arma;

class SGD {
public:
    SGD(System* system, double eta);
    void SGDOptimize(arma::vec gradE);

private:
    double m_eta;   // learning rate
    arma::vec a;    // visible bias
    arma::vec b;    // hidden bias
    arma::mat w;    // interaction matrix between biases

    class System* m_system = nullptr;

};