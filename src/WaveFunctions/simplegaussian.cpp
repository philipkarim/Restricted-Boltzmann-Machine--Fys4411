#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>
#include <armadillo>

using namespace std;
SimpleGaussian::SimpleGaussian(System* system, int part, 
                                int dim, double sigma) :
        WaveFunction(system) {
    m_sigma = sigma;
    m_system->setNumberParticles(part);
    m_system->setNumberDimensions(dim);

}

double SimpleGaussian::evaluate(arma::vec position) {
     //Implementation of wavefunction at the given position


    //Return a double value
    return 1.;
}

double SimpleGaussian::computeDoubleDerivative(arma::vec position) {
     //Computes the value of the analytical double derivative for the non interacting case. 

    //Return a double value
    return 1.;
  }

double SimpleGaussian::computeDerivative(arma::vec position) {
     //Computes the value of the analytical derivative 

    //Return a double value
    return 1.0;
  }