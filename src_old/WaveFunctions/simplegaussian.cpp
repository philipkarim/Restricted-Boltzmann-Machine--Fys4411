#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>
#include <vector>

using namespace std;
SimpleGaussian::SimpleGaussian(System* system, int n_hidden, int n_visible, int part, int dim, double sigma, bool gaussian, double initialization) :
        WaveFunction(system) {
    assert(n_hidden >= 0);
    assert(n_visible >= 0);
    m_nh = n_hidden;
    m_nv = n_visible;
    m_sigma = sigma;
    m_gaussian = gaussian;
    m_initialization = initialization;

    m_system->setNumberHiddenNodes(n_hidden);
    m_system->setNumberVisibleNodes(n_visible);
    m_system->setNumberParticles(part);
    m_system->setNumberDimensions(dim);

    //Initialize 
    setupInitialState();

}

double SimpleGaussian::evaluate(std::vector<double> position) {
     //Implementation of wavefunction at the given position


    //Return a double value
    return 1.
}

void SimpleGaussian::setupInitialState(){
    m_x.zeros(m_nv);
    m_h.zeros(m_nh);
    m_a.zeros(m_nv);
    m_b.zeros(m_nh);
    m_w.zeros(m_nv, m_nh);

    std::uniform_real_distribution<double> uniform_weights(-m_initialization, m_initialization);
    std::uniform_real_distribution<double> uniform_position(-0.5, 0.5);
    std::normal_distribution<double> normal_weights(0, m_initialization);

    // initialize weights according to either a uniform or gaussian distribution
    if (m_gaussian == "Uniform"){
        for (int i=0; i<m_nv; i++){
            m_a[i] = uniform_weights(m_randomEngine);
            for (int j=0; j<m_nh; j++){
                m_w(i, j) = uniform_weights(m_randomEngine);
            }
        }

        for (int j=0; j<m_nh; j++){
            m_b[j] = uniform_weights(m_randomEngine);
        }
    }

    else if (m_gaussian == "Normal"){
        for (int i=0; i<m_nv; i++){
            m_a[i] = normal_weights(m_randomEngine);
            for (int j=0; j<m_nh; j++){
                m_w(i, j) = normal_weights(m_randomEngine);
            }
        }

        for (int j=0; j<m_nh; j++){
            m_b[j] = normal_weights(m_randomEngine);
        }
    }

    for (int i=0; i<m_nv; i++){
        m_x[i] = uniform_position(m_randomEngine);
    }

}

double SimpleGaussian::computeDoubleDerivative(std::vector<double> position) {
     //Computes the value of the analytical double derivative for the non interacting case. 

    //Return a double value
    return 1.
  }

double SimpleGaussian::computeDerivative(std::vector<double> position) {
     //Computes the value of the analytical derivative 

    //Return a double value
    return 1.0
  }