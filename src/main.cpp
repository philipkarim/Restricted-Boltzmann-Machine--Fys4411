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

    int numberOfSteps       = (int) pow(2,20); //Amount of metropolis steps
    int cycles_RBM          = 100;
    int numberOfDimensions  = 1;            // Set amount of dimensions
    int numberOfParticles   = 1;            // Set amount of particles
    int hidden_nodes        = 2;
    int visible_nodes       = numberOfDimensions*numberOfParticles;
    int sampler_method      = 0;            //0=BF, 1=IS, 2=GS
    bool uniform_distr      = false;         //Normal=false, Uniform=true
    double omega            = 1.0;          // Oscillator frequency.
    double stepLength       = 0.5;          // Metropolis step length.
    double timeStep         = 0.25;         // Metropolis time step (Importance sampling)
    double equilibration    = 0.2;          // Amount of the total steps used for equilibration.
    bool interaction        = false;        // True-> interaction, False->Not interaction
    double sigma_val        = 1.;
    double initialization   = 0.001;
    double learningRate     = 0.001;
    //Write to file
    bool generalwtf        =false;          // General information- write to file
    bool explore_distribution=true;
    bool find_optimal_step =false;

    System* system = new System(seed);
    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setWaveFunction             (new NeuralState(system, numberOfParticles, numberOfDimensions, sigma_val));
    system->setInitialState             (new RandomUniform(system, hidden_nodes, visible_nodes, uniform_distr, initialization));
    system->setStepLength               (stepLength);                                                               
    system->setLearningRate             (learningRate);
    system->setTimeStep                 (timeStep);
    system->setEquilibrationFraction    (equilibration);
    system->setSampleMethod             (sampler_method);
    system->setInteraction              (interaction);
    system->setgeneralwtf               (generalwtf);
    
    if(explore_distribution==true){
      int pid, pid1, pid2, pid3, pid4, pid5, pid6;
      //Using more cores to achieve more results faster
      pid=fork();
      if(pid==0){
        //Uniform distribution, initialization value=0.25
          system->setHamiltonian              (new HarmonicOscillator(system, omega));
          system->setWaveFunction             (new NeuralState(system, numberOfParticles, numberOfDimensions, sigma_val));
          system->setInitialState             (new RandomUniform(system, hidden_nodes, visible_nodes,true, 0.25)); 
          system->setWtfDistibution           (explore_distribution);
          system->runBoltzmannMachine         (cycles_RBM, numberOfSteps);}

        //Uniform distribution, initialization value=0.01
      else{pid1=fork(); if(pid1==0){
          system->setHamiltonian              (new HarmonicOscillator(system, omega));
          system->setWaveFunction             (new NeuralState(system, numberOfParticles, numberOfDimensions, sigma_val));
          system->setInitialState             (new RandomUniform(system, hidden_nodes, visible_nodes,true, 0.01)); 
          system->setWtfDistibution           (explore_distribution);
          system->runBoltzmannMachine         (cycles_RBM, numberOfSteps);}

        //Uniform distribution, initialization value=0.001
      else{pid2=fork(); if(pid2==0){
          system->setHamiltonian              (new HarmonicOscillator(system, omega));
          system->setWaveFunction             (new NeuralState(system, numberOfParticles, numberOfDimensions, sigma_val));
          system->setInitialState             (new RandomUniform(system, hidden_nodes, visible_nodes,true, 0.001)); 
          system->setWtfDistibution           (explore_distribution);
          system->runBoltzmannMachine         (cycles_RBM, numberOfSteps);}

        //Normal distribution, initialization value=0.001
      else{pid3=fork(); if(pid3==0){
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new NeuralState(system, numberOfParticles, numberOfDimensions, sigma_val));
        system->setInitialState             (new RandomUniform(system, hidden_nodes, visible_nodes,true, 0.005)); 
        system->setWtfDistibution           (explore_distribution);
        system->runBoltzmannMachine         (cycles_RBM, numberOfSteps);}

        //Normal distribution, initialization value=0.001
      else{pid4=fork(); if(pid4==0){
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new NeuralState(system, numberOfParticles, numberOfDimensions, sigma_val));
        system->setInitialState             (new RandomUniform(system, hidden_nodes, visible_nodes,false, 0.25)); 
        system->setWtfDistibution           (explore_distribution);
        system->runBoltzmannMachine         (cycles_RBM, numberOfSteps);}

      else{pid5=fork(); if(pid5==0){
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new NeuralState(system, numberOfParticles, numberOfDimensions, sigma_val));
        system->setInitialState             (new RandomUniform(system, hidden_nodes, visible_nodes,false, 0.01)); 
        system->setWtfDistibution           (explore_distribution);
        system->runBoltzmannMachine         (cycles_RBM, numberOfSteps);}

      else{pid6=fork(); if(pid6==0){
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new NeuralState(system, numberOfParticles, numberOfDimensions, sigma_val));
        system->setInitialState             (new RandomUniform(system, hidden_nodes, visible_nodes,false, 0.001)); 
        system->setWtfDistibution           (explore_distribution);
        system->runBoltzmannMachine         (cycles_RBM, numberOfSteps);}

        //Normal distribution, initialization value=0.001
      else{
          system->setHamiltonian              (new HarmonicOscillator(system, omega));
          system->setWaveFunction             (new NeuralState(system, numberOfParticles, numberOfDimensions, sigma_val));
          system->setInitialState             (new RandomUniform(system, hidden_nodes, visible_nodes,false, 0.005)); 
          system->setWtfDistibution           (explore_distribution);
          system->runBoltzmannMachine         (cycles_RBM, numberOfSteps);
        }}}}}}}}

    else if(find_optimal_step==true){
      int pid, pid1, pid2;
      //Using more cores to achieve more results faster
      pid=fork();
      if(pid==0){
        //bf, non interacting, different step sizes
        for (double i=1.5; i>0.1; i-=0.1){
          system->setHamiltonian              (new HarmonicOscillator(system, omega));
          system->setWaveFunction             (new NeuralState(system, 2, 2, sigma_val));
          system->setInitialState             (new RandomUniform(system, 2, 4,uniform_distr, initialization)); 
          system->setStepLength               (i);                                                   
          system->setLearningRate             (learningRate);
          system->setTimeStep                 (timeStep);
          system->setEquilibrationFraction    (equilibration);
          system->setSampleMethod             (0);
          system->setInteraction              (true);
          system->setgeneralwtf               (generalwtf);
          system->setwtfSteps                 (find_optimal_step);
          system->runBoltzmannMachine         (cycles_RBM, numberOfSteps);
        }}

      /*
      //bf, interacting, different step sizes
      else{pid1=fork(); if(pid1==0){
        for (double k=1.5; k>0.1; k-=0.1){
          system->setHamiltonian              (new HarmonicOscillator(system, omega));
          system->setWaveFunction             (new NeuralState(system, 2, 2, sigma_val));
          system->setInitialState             (new RandomUniform(system, 2, 4,uniform_distr, initialization));
          system->setStepLength               (k);                                                        
          system->setLearningRate             (learningRate);
          system->setTimeStep                 (timeStep);
          system->setEquilibrationFraction    (equilibration);
          system->setSampleMethod             (0);
          system->setInteraction              (true);
          system->setgeneralwtf               (generalwtf);
          system->setwtfSteps                 (find_optimal_step);
          system->runBoltzmannMachine         (cycles_RBM, numberOfSteps);
        }}

      //is, non interacting, different timesteps
      else{pid2=fork(); if(pid2==0){
        for (double j=1; j>0.1; j-=0.1){
          system->setHamiltonian              (new HarmonicOscillator(system, omega));
          system->setWaveFunction             (new NeuralState(system, 1, 1, sigma_val));
          system->setInitialState             (new RandomUniform(system, 2, 1, uniform_distr, initialization));
          system->setStepLength               (stepLength);                                                  
          system->setLearningRate             (learningRate);
          system->setTimeStep                 (j);
          system->setEquilibrationFraction    (equilibration);
          system->setSampleMethod             (1);
          system->setInteraction              (false);
          system->setgeneralwtf               (generalwtf);
          system->setwtfSteps                 (find_optimal_step);
          system->runBoltzmannMachine         (cycles_RBM, numberOfSteps);
        }}
        */
      else{
        for (double l=1; l>0.1; l-=0.1){
          system->setHamiltonian              (new HarmonicOscillator(system, omega));
          system->setWaveFunction             (new NeuralState(system, 2, 2, sigma_val));
          system->setInitialState             (new RandomUniform(system, 2, 4, uniform_distr, initialization));
          system->setStepLength               (stepLength);                                                     
          system->setLearningRate             (learningRate);
          system->setTimeStep                 (l);
          system->setEquilibrationFraction    (equilibration);
          system->setSampleMethod             (1);
          system->setInteraction              (true);
          system->setgeneralwtf               (generalwtf);
          system->setwtfSteps                 (find_optimal_step);
          system->runBoltzmannMachine         (cycles_RBM, numberOfSteps);
        }}} //}}
    else{
    //Setting the different values defined higher in the code
    system->runBoltzmannMachine         (cycles_RBM, numberOfSteps);
    }
    return 0;
}

//Tasks:
//-----------------------
//See if it handles maximum of steps and time steps

//If all results makes sense it is time to extract the results, compute the results needed:
  //Write to files the things we want
  //Compute the energy by blocking method
  //Plot things



/*
For compilers to find openblas you may need to set:
  export LDFLAGS="-L/usr/local/opt/openblas/lib"
  export CPPFLAGS="-I/usr/local/opt/openblas/include"
*/