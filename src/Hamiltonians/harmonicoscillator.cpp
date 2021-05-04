#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

#include <vector>

using std::cout;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega   > 0);
    m_omega  = omega;
}

double HarmonicOscillator::computeLocalEnergy(std::vector<double> X_visible) {
    //This function is computing the kinetic and potential energies

    //Defining some variables to be used in the calculations 
    double potentialEnergy = 0;
    double kineticEnergy   = 0;
    double interactionEnergy=0;
    //Evaluating the wave function

    //Computing the non interacting energy
    if (m_system->getInteraction()==false){

    //Computing the kinetic energy
    //kineticEnergy=m_system->getWaveFunction()->computeDoubleDerivative();

    potentialEnergy=computePotentialEnergy(X_visible);

    //Returning the energy as a double
    return (kineticEnergy + potentialEnergy);
  }
    else{
      cout<<"Interacting energy not done";
    
    return 1.;
    }
  
  }
   
double HarmonicOscillator::computePotentialEnergy(std::vector<double> X_visible) {
  //Done

  //Potential energy noninteracting case
  int numberVN=0; //= m_system->getNumberVisibleNodes();
  //Defining some variables to be used later
  double potentialEnergy=0;
    
    //    int particleCounter = 0;
    for(int i = 0; i<numberVN; i++){
        potentialEnergy += X_visible[i]*X_visible[i];
    }

    //Returning the potential energy
    return potentialEnergy*0.5*m_omega*m_omega;

}

