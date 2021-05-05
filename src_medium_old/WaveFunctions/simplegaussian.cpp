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
    return 1.;
}

void SimpleGaussian::setupInitialState(){

}

double SimpleGaussian::computeDoubleDerivative(std::vector<double> position) {
     //Computes the value of the analytical double derivative for the non interacting case. 

    //Return a double value
    return 1.;
  }

double SimpleGaussian::computeDerivative(std::vector<double> position) {
     //Computes the value of the analytical derivative 

    //Return a double value
    return 1.0;
  }