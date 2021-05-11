#pragma once
#include "hamiltonian.h"
#include <vector>
#include <armadillo>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);
    double computeLocalEnergy(arma::vec X_visible);
    double computePotentialEnergy(arma::vec X_visible);
    double computeKineticEnergy(arma:: vec X_visible);
    double computeInteractingEnergy(arma:: vec X_visible);
    double particleDistance(int i, int j, arma:: vec X_visible);



private:
    double m_omega = 0;


};
