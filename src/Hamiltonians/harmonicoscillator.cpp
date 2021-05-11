#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

#include <armadillo>

using std::cout;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
    Hamiltonian(system) {
    assert(omega   > 0);
    m_omega  = omega;
    
}

double HarmonicOscillator::computeLocalEnergy(arma::vec X_visible) {
  //This function is computing the kinetic and potential energies

  //Defining some variables to be used in the calculations 
  double potentialEnergy = 0;
  double kineticEnergy   = 0;
  double interactionEnergy=0;
  //Evaluating the wave function

  kineticEnergy=computeKineticEnergy(X_visible);
  potentialEnergy=computePotentialEnergy(X_visible);

  if (m_system->getInteraction()==true){
    interactionEnergy=computeInteractingEnergy(X_visible);
  }
  else{
    interactionEnergy=0;
    }
  
  return kineticEnergy+potentialEnergy+interactionEnergy;
  }
   
double HarmonicOscillator::computePotentialEnergy(arma::vec X_visible) {
  //Potential energy
  //Defining some variables to be used
  double potentialEnergy=0;
    
    //    int particleCounter = 0;
    for(int i = 0; i<m_system->getNumberOfVN(); i++){
        potentialEnergy += X_visible[i]*X_visible[i];
    }

    //Returning the potential energy
    return potentialEnergy*0.5*m_omega*m_omega;

}

double HarmonicOscillator::computeKineticEnergy(arma:: vec X_visible){
  double double_derivative, derivative;
  double_derivative=m_system->getWaveFunction()->computeDoubleDerivative(X_visible);
  derivative       =m_system->getWaveFunction()->computeDerivative(X_visible);

  return -0.5*(derivative*derivative+double_derivative);

}

double HarmonicOscillator::computeInteractingEnergy(arma:: vec X_visible){
  int n_particles=m_system->getNumberOfParticles();
  double distance;
  double inter_energ;
  for (int i_1=0; i_1<n_particles-1; i_1++){
    for(int i_2=i_1+1; i_2<n_particles; i_2++){
      distance = particleDistance(i_1,i_2, X_visible);
      inter_energ+= 1./distance;
    }
  }  
  return inter_energ;

}

double HarmonicOscillator::particleDistance(int i, int j, arma:: vec X_visible){
    int dims = m_system->getNumberOfDimensions();
    double distance2=0;

    for(int d = 0; d<dims; d++){
        distance2+= (X_visible[dims*i+d]-X_visible[dims*j+d])*(X_visible[dims*i+d]-X_visible[dims*j+d]);
    }
    return sqrt(distance2);
}

//Harmonic Oscillator is done!