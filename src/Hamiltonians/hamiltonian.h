#pragma once
#include <armadillo>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(arma::vec X_visible)=0;
    virtual double computePotentialEnergy(arma::vec X_visible)=0;

protected:
    class System* m_system = nullptr;
};
