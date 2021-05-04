#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include <iostream>

#include "WaveFunctions/simplegaussian.h"

#include <stdlib.h>     /* exit, EXIT_FAILURE */


using namespace std;

System::System() {
    m_random = new Random();
}

System::System(int seed) {
    m_random = new Random(seed);
}

bool System::metropolisStep() {
    // Performing the actual Metropolis step for the Metropolis algorithm:
    // Choosing a particle at random and changing it's position by a random
    // amount, and checks if the step is accepted by the Metropolis test
    
    // Defining some variables to be used
     int random_index;
     double psi_factor, step;
     double wfold=m_waveFunction->evaluate(m_particles);
     std::vector<double> PositionOld=std::vector<double>();

     //Random integer generator
     std::random_device rd;
     std::mt19937_64 gen(rd());
     std::uniform_int_distribution<int> distribution(0,m_numberOfParticles-1);
     std::uniform_real_distribution<double> UniformNumberGenerator(0.0,1.0);

     //Random index used to choose a random particle
     random_index=distribution(gen);
     //Defining the random particle:
     PositionOld=m_particles[random_index]->getPosition();

     //Start the step which gives movement to the particle
     for (int dim=0; dim<m_numberOfDimensions; dim++){
        step=m_stepLength*(UniformNumberGenerator(gen)-0.5);
        m_particles[random_index]->adjustPosition(step, dim);
     }

     //Extracting the new wavefunction, and checks if it is accepted
     double wfnew=m_waveFunction->evaluate(m_particles);
     psi_factor=wfnew*wfnew/(wfold*wfold);
     
     //Checks if the move is accepted:
     if (UniformNumberGenerator(gen)<=psi_factor){
        wfold=wfnew;
        return true;
     }
     else{
         m_particles[random_index]->setPosition(PositionOld);
        return false;
      }
}

bool System::metropolisStepImportanceSampling() {
    // Performing the actual Metropolis step for the Metropolis- Hastings
    //algorithm: Choosing a particle at random and changing it's position 
    // by a random amount, and checks if the step is accepted by the 
    //Metropolis-Hastings test
    
    return true

}

bool System::GibbsSampling() {
    // Performing the actual Metropolis step for the Metropolis- Hastings
    //algorithm: Choosing a particle at random and changing it's position 
    // by a random amount, and checks if the step is accepted by the 
    //Metropolis-Hastings test
    
    return true

}

void System::runMetropolisSteps(int RBMCycles, int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_RBMCycles                 = RBMCycles;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    bool acceptedStep;
    
    //Looping over the amount of metropolis steps
    //for either the Metroopolis algorithm or the
    //Metropolis-Hastings algorithm
    for (int i=0; i < numberOfMetropolisSteps; i++) {
      if (m_sampleMethod=="Brute Force"){
        acceptedStep = metropolisStep();}
      else if(m_sampleMethod=="Importance Sampling"){
        acceptedStep = metropolisStepImportanceSampling();
      }
      else if(m_sampleMethod=="Gibbs Sampling"){
        acceptedStep = metropolisStepImportanceSampling();
      }
      else{
          cout<<"---No sampling method chosen---";
          exit (EXIT_FAILURE);
      }

      //If statement to send the accepted steps into the sampler
      //after the system is at rest
      if (i>=numberOfMetropolisSteps*m_equilibrationFraction){
        m_sampler->sample(acceptedStep);
      }
    //cout<<"The current step is "<<i<<endl;
    }

    //Chooosing what to sample
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
    if (m_general_wtf==true){m_sampler->writeToFile();}
    
}



void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setSampleMethod(string sampleMethod) {
    m_sampleMethod = sampleMethod;
}

void System::setTimeStep(double timeStep) {
    m_timeStep= timeStep;
}

void System::setInteraction(bool interaction) {
    m_interaction= interaction;
}

void System::setgeneralwtf(bool generalwtf) {
    m_general_wtf = generalwtf;
}

void System::setNumberParticles(int numberOfParticles) {
    m_numberParticles = numberOfParticles;
}

void System::setNumberDimensions(int numberOfDimensions) {
    m_numberDimensions = numberOfDimensions;
}
