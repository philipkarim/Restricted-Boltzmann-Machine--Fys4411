#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);
    double computeLocalEnergy(arma::vec X_visible);
    double computePotentialEnergy(arma::vec X_visible);
private:
    double m_omega = 0;


};
