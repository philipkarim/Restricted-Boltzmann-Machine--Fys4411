#include "SGD.h"
#include <cassert>
#include <iostream>
#include <armadillo>
#include "system.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"


SGD::SGD(System* system, double eta){
    m_eta = eta;
    m_system = system;

    m_system->setLearningRate(m_eta);
}

/* Calculate optimized weights via stochastic gradient descent */
void SGD::SGDOptimize(arma::vec gradE){
    int m_nv = m_system->getNumberOfVN();
    int m_nh = m_system->getNumberOfHN();

    // get current biases
    a = m_system->getWaveFunction()->get_a();
    b = m_system->getWaveFunction()->get_b();
    w = m_system->getWaveFunction()->get_w();

    // optimize biases
    for (int i=0; i<m_nv; i++){
        a[i] = a[i] - m_eta*gradE[i];
    }

    for (int j=0; j<m_nh; j++){
        b[j] = b[j] - m_eta*gradE[m_nv + j];

    }

    int count = m_nh + m_nv;
    for (int i=0; i<m_nv; i++){
        for (int j=0; j<m_nh; j++){
            w(i, j) = w(i, j) - m_eta*gradE[count];
            count++;
        }
    }

    // set optimized biases as new weights
    m_system->getWaveFunction()->set_a(a);
    m_system->getWaveFunction()->set_b(b);
    m_system->getWaveFunction()->set_w(w);
}