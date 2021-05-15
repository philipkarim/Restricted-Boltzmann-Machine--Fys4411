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
using namespace arma;

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

    //cout<<step;
    
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

    //Declaring vaiables to be used:
    double Position_old, wfnew, wfold, part_1, part_2, rand_norm ,green_factor, step, greenRate=0;
    int random_index;
    //Defining position and quantum force vectors
    //to be used in the importance sampling
    vec X_old=m_waveFunction->get_X();
    double QFOld;
    vec X_new=m_waveFunction->get_X();
    double QFNew;

    //Random integer generator
    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<int> distribution(0,m_numberOfVN-1);//-1?
    uniform_real_distribution<double> UniformNumberGenerator(0.0,1.0);
    normal_distribution<double> Normaldistribution(0.0,1.0);

    random_index  =distribution(gen);
    rand_norm  =Normaldistribution(gen);
    
    //Defining the values of the previous position
    wfold=m_waveFunction->evaluate(X_old);    
    QFOld=m_waveFunction->computeQuantumForce(X_old, random_index);
    Position_old= X_new[random_index];
    X_new[random_index]+=QFOld*m_timeStep*0.5 + sqrt(m_timeStep)*rand_norm;

    // Evaluate new quantities
    wfnew = m_waveFunction->evaluate(X_new);
    QFNew=m_waveFunction->computeQuantumForce(X_new, random_index);

    // Compute greens function
    part_1=X_old[random_index]-X_new[random_index]-0.5*m_timeStep*QFNew;
    part_2=X_new[random_index]-X_old[random_index]-0.5*m_timeStep*QFOld;
    greenRate=(part_2*part_2)-(part_1*part_1);

    greenRate = exp(greenRate/(2*m_timeStep));
    green_factor = greenRate*wfnew*wfnew/(wfold*wfold);

    // Check if the step is accepted
    if (UniformNumberGenerator(gen) <= green_factor) {
        m_waveFunction->set_X(X_new);
        return true;
    }
    else {
        //Print the random index here and over, to see if it is the same index
        X_old[random_index]=Position_old;
        return false;
    }

}

bool System::GibbsSampling() {
    // Performing Gibbs sampling
    
    double randu, P, sigma;
   vec O;

    sigma = m_waveFunction->getSigma();
    vec X = m_waveFunction->get_X();
    vec m_h = m_waveFunction->get_h();
    vec m_a = m_waveFunction->get_a();
    vec m_b = m_waveFunction->get_b();
    mat m_w = m_waveFunction->get_w();

    // Get probability, P(h=1|x), of hidden values to equal 1, given the logistic sigmoid function
    // Set hidden values equal to 0 if probability less than a randuom uniform variable
    O = m_b + ((X.t()*m_w).t()) / (sigma*sigma);
    for (int j=0; j < m_numberHiddenNodes; j++){
        randu = getUniform(0, 1);
        P = 1.0/(1+exp(-O[j]));     // probability from logistic sigmoid
        if (P < randu){ m_h[j] = 0; }
        else{ m_h[j] = 1; }
    }

    // Set the new positions according to the hidden nodes
    double randn, x_mean;
    vec new_pos, wh;
    new_pos.zeros(m_numberVisibleNodes);

    wh = m_w*m_h;
    for (int i=0; i<m_numberVisibleNodes; i++){
        x_mean = m_a[i] + wh[i];
        new_pos[i] = getGaussian(x_mean, sigma);
    }

//    cout << "X before: " << endl; m_waveFunction->get_X().print();
    m_waveFunction->set_X(new_pos);
//    cout << "X after: " << endl; m_waveFunction->get_X().print();

    return true;
}
    

void System::runBoltzmannMachine(int RBMCycles, int numberOfMetropolisSteps, double lr){
    m_RBMCycles                 = RBMCycles;
    m_SGD= (new SGD(this, lr));

    //m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    //m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (int rbm_cycle=0; rbm_cycle<RBMCycles; rbm_cycle++){
        cout<<"RBM cycle: "<<rbm_cycle<<"\n------------------------------";
        m_sampler                   = new Sampler(this);
        //m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
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
