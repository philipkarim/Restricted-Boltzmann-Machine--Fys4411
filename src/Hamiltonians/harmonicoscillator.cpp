#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

#include <armadillo>

using namespace std;
using namespace arma;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
    Hamiltonian(system) {
    assert(omega   > 0);
    m_omega  = omega;
    
}

double HarmonicOscillator::computeLocalEnergy() {
  //This function is computing the kinetic and potential energies
  vec X_visible=m_system->getWaveFunction()->get_X();

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
   
double HarmonicOscillator::computePotentialEnergy(vec X_visible) {
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

double HarmonicOscillator::computeKineticEnergy(vec X_visible){
  double double_derivative, derivative;
  if(m_system->getSampleMethod()==2){
    double_derivative=0.5*m_system->getWaveFunction()->computeDoubleDerivative(X_visible);
    derivative       =0.5*m_system->getWaveFunction()->computeDerivative(X_visible);
  }
  else{
    double_derivative=m_system->getWaveFunction()->computeDoubleDerivative(X_visible);
    derivative       =m_system->getWaveFunction()->computeDerivative(X_visible);
  }

  return -0.5*(derivative*derivative+double_derivative);

}

double HarmonicOscillator::computeInteractingEnergy( vec X_visible){
  //Computing the interacting energy, sends in two particles
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

double HarmonicOscillator::particleDistance(int i, int j,  vec X_visible){
    int dims = m_system->getNumberOfDimensions();
    double norm=0;
    double part_1, part_2;
    for(int d = 0; d<dims; d++){
        part_1=X_visible[dims*i+d]-X_visible[dims*j+d];
        part_2=X_visible[dims*i+d]-X_visible[dims*j+d];
        norm+= part_1*part_2;
    }
    return sqrt(norm);
}

//Harmonic Oscillator is done down to here:

//Computes the derivative with respect to the different parameters, to be used in SGD
vec HarmonicOscillator::computeParameterDerivatives(){
    //Defining some variables
    vec xx = m_system->getWaveFunction()->get_X();
    vec aa = m_system->getWaveFunction()->get_a();
    vec bb = m_system->getWaveFunction()->get_b();
    mat ww = m_system->getWaveFunction()->get_w();
    double sig = m_system->getWaveFunction()->getSigma();

    int numberOfVN = m_system->getNumberOfVN();
    int numberOfHN = m_system->getNumberOfHN();

    //Continue under and change name of the function to something easier
    //Double check the formulas to get rid of uneccesarry computations .t
    vec O = bb + (xx.t()*ww).t()*(1/((double) sig*sig));
    vec dPsi; dPsi.zeros(numberOfVN + numberOfHN + numberOfVN*numberOfHN);

    // compute d(psi)/d(a)*1/psi
    for (int i=0; i<numberOfVN; i++){
        dPsi[i] = (xx[i] - aa[i])/(sig*sig);
    }

    // compute d(psi)/d(b)*1/psi
    for (int i=numberOfVN; i<numberOfVN+numberOfHN; i++){
        dPsi[i] = 1.0/(exp(-O[i - numberOfVN]) + 1);
    }

    // compute d(psi)/d(w)*1/psi
    int i = numberOfVN + numberOfHN;
    for (int j=0; j<numberOfVN; j++){
        for (int k=0; k<numberOfHN; k++){
            dPsi[i] = xx[j]/(exp(-O[k]) + 1)*(1.0/sig*sig);
            i++;
        }
    }

    return dPsi;
}
