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
    bool uniform_distr      = 0;            //0=Normal, 1=Uniform
    double omega            = 1.0;          // Oscillator frequency.
    double stepLength       = 0.5;          // Metropolis step length.
    double timeStep         = 0.25;         // Metropolis time step (Importance sampling)
    double equilibration    = 0.2;          // Amount of the total steps used for equilibration.
    bool interaction        = false;        // True-> interaction, False->Not interaction
    double sigma_val        =1.0;
    double initialization   =0.001;
    double learningRate     =0.001;
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
    system->runBoltzmannMachine          (cycles_RBM, numberOfSteps, learningRate);

    return 0;
}

//Tasks:
//-----------------------
//Read through all code
//Implement importance and gibbs and gibbs energy