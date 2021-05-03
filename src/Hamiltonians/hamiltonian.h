#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(std::vector<double> X_visible)=0;
    virtual double computePotentialEnergy(std::vector<double> X_visible)=0;

protected:
    class System* m_system = nullptr;
};
