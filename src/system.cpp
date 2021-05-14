#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include <iostream>
#include "SGD.h"
#include<armadillo>
#include "WaveFunctions/neuralstate.h"

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

    //Random integer generator
    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<int> distribution(0,m_numberOfVN-1);//-1?
    uniform_real_distribution<double> UniformNumberGenerator(0.0,1.0);

    // Defining some variables to be used
    int random_index;
    double Position_old, psi_old, psi_new, psi_factor, step;
    arma::vec X_old;
    
    X_old=m_waveFunction->get_X();

    //Random index used to choose a random particle
    random_index=distribution(gen);

    //Defining the random particle:
    Position_old=X_old[random_index];
    psi_old=m_waveFunction->evaluate(X_old);

    //Start the step which gives movement to the particle
    step=m_stepLength*(UniformNumberGenerator(gen)-0.5);
    X_old[random_index]+=step;
    
     //Extracting the new wavefunction, and checks if it is accepted
    psi_new=m_waveFunction->evaluate(X_old);
    psi_factor=psi_new*psi_new/(psi_old*psi_old);
     
    //Checks if the move is accepted:
    if (UniformNumberGenerator(gen)<=psi_factor){
        m_waveFunction->set_X(X_old);        

        return true;
     }
     else{
         X_old[random_index]=Position_old;
        
        return false;
     }
}

bool System::metropolisStepImportanceSampling() {
    // Performing the actual Metropolis step for the Metropolis- Hastings
    //algorithm: Choosing a particle at random and changing it's position 
    // by a random amount, and checks if the step is accepted by the 
    //Metropolis-Hastings test
    
    return true;

}

bool System::GibbsSampling() {
    // Performing the actual Metropolis step for the Metropolis- Hastings
    //algorithm: Choosing a particle at random and changing it's position 
    // by a random amount, and checks if the step is accepted by the 
    //Metropolis-Hastings test
    
    return true;

}

void System::runBoltzmannMachine(int RBMCycles, int numberOfMetropolisSteps, double lr){
    m_RBMCycles                 = RBMCycles;
    m_SGD= (new SGD(this, lr));

    //m_sampler                   = new Sampler(this);
    //m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    //m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (int rbm_cycle=0; rbm_cycle<RBMCycles; rbm_cycle++){
        m_sampler                   = new Sampler(this);
        m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
        m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
        runMetropolisSteps();
        m_SGD->SGDOptimize(m_sampler->getGradient());
    }


}
void System::runMetropolisSteps() {
    bool acceptedStep;
    
    //Looping over the amount of metropolis steps
    //for either of the sampling methods
    for (int i=0; i < m_numberOfMetropolisSteps; i++) {
      if (m_sampleMethod==0){
        acceptedStep = metropolisStep();}
      else if(m_sampleMethod==1){
        acceptedStep = metropolisStepImportanceSampling();
      }
      else if(m_sampleMethod==2){
        acceptedStep = metropolisStepImportanceSampling();
      }
      else{
          cout<<"---No sampling method chosen---";
          exit (EXIT_FAILURE);
      }

      //If statement to send the accepted steps into the sampler
      //after the system is at rest
      if (i>=m_numberOfMetropolisSteps*m_equilibrationFraction){
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

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}
/*
void System::setOptimizer(SGD* learningRate) {
    m_learningRate = learningRate;
}
*/
void System::setSampleMethod(int sampleMethod) {
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
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setLearningRate(double learningRate) {
    m_learningRate = learningRate;
}

void System::setNumberOfHN(int n_hidden){ 
    m_numberOfHN = n_hidden; 
}

void System::setNumberOfVN(int n_visible){ 
    m_numberOfVN = n_visible; 
}
