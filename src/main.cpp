#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/neuralstate.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "SGD.h"
#include <unistd.h>

#include <string>

using namespace std;

int main() {

    // Seed for the random number generator
    int seed = 2021;

    //Correct approx: 12, 10, 1, 1, 4, 0.03
    //Correct approx: 18, 50, 1, 1, 4, 0.001


    int numberOfSteps       = (int) pow(2,13); //Amount of metropolis steps
    int cycles_RBM          = 20;
    int numberOfDimensions  = 2;            // Set amount of dimensions
    int numberOfParticles   = 2;            // Set amount of particles
    int hidden_nodes        = 2;
    int visible_nodes       = numberOfDimensions*numberOfParticles;
    int sampler_method      = 2;            //0=BF, 1=IS, 2=GS
    bool uniform_distr      = 1;//Is normal only for gibbs?            //0=Normal, 1=Uniform
    double omega            = 1.0;          // Oscillator frequency.
    double stepLength       = 0.5;          // Metropolis step length.
    double timeStep         = 0.25;         // Metropolis time step (Importance sampling)
    double equilibration    = 0.2;          // Amount of the total steps used for equilibration.
    bool interaction        = false;        // True-> interaction, False->Not interaction
    double sigma_val        = 1.0;
    double initialization   = 0.001;
    double learningRate     = 0.001;
    //Write to file
    bool generalwtf        =false;          // General information- write to file

    //Setting the different values defined higher in the code
    System* system = new System(seed);
    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setWaveFunction             (new NeuralState(system, numberOfParticles, numberOfDimensions, sigma_val));
    system->setInitialState             (new RandomUniform(system, hidden_nodes, visible_nodes, 
                                                                   uniform_distr, initialization));
    system->setLearningRate             (learningRate);
    system->setStepLength               (stepLength);
    system->setTimeStep                 (timeStep);
    system->setEquilibrationFraction    (equilibration);
    system->setSampleMethod             (sampler_method);
    system->setInteraction              (interaction);
    system->setgeneralwtf               (generalwtf);
    system->runBoltzmannMachine         (cycles_RBM, numberOfSteps);

    return 0;
}

//Tasks:
//-----------------------
//Read through all code-->Done almost
//Implement gibbs sampling-->Done
//Implement Gibbs energy--> Done but not sure if it is correct
//Move sigmoid in neural.h from private to public?
//Update quantum force the same way is in derivative since both is actually transposed
//Try and see if hidden nodes are updates now with different distribution

//If all results makes sense it is time to extract the results, compute the results needed:
  //Write to files the things we want
  //Compute the energy by blocking method
  //Plot things
/*
Main—>Done
system—>Evereything done but Gibbs sampling 
Sampler—> Done
Sgd—>Done-Just check a small thing
Hamiltonians/—> everything done but the last function, need to check the pdf to get rid of transposes

Wavefunctions/-> Done but QForce
Initial states/—> Done
*/


/*
For compilers to find openblas you may need to set:
  export LDFLAGS="-L/usr/local/opt/openblas/lib"
  export CPPFLAGS="-I/usr/local/opt/openblas/include"
*/