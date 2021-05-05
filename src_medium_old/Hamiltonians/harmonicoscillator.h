#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);
    double computeLocalEnergy(std::vector<double> X_visible);
    double computePotentialEnergy(std::vector<double> X_visible);
private:
    double m_omega = 0;


};
