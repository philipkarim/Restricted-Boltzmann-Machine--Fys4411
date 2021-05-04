#include <iostream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
 #include <iomanip>

//Write to file modules
#include <fstream>
//Modules to meassure the cpu time
#include <ctime>
#include <ratio>
#include <chrono>

#include <string>

using namespace std::chrono;
using namespace std;

using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep) {
    //Sampling interesting results
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_cumulativeEnergy2 = 0;
        time_sec =0;
        m_cumulativeE_Lderiv=0;
        m_cumulativeE_Lderiv_expect=0;
    }


    
    //Starting the clock
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    //Calculating the local energy
    double localEnergy; //= m_system->getHamiltonian()->computeLocalEnergy(hidden_nodes_maybe);

   //Stopping the clock, adding the time together for each energy cycle
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
	  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    time_sec += time_span.count();

    //Saving values to be used in blocking
    //Metropolis steps are 2^n, performs blocking on 2^(n-1) because of the equilibration factor
    if (meanenergy_list.size()<m_system->getNumberOfMetropolisSteps()/2){
      meanenergy_list.push_back(localEnergy);
    }

    //Cumulating the energy
    m_cumulativeEnergy  += localEnergy;
    m_stepNumber++;

    //Used in variance
    m_cumulativeEnergy2+=(localEnergy*localEnergy);

    if (acceptedStep){
        m_acceptedSteps++;
    }

}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    double  ef = m_system->getEquilibrationFraction();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " CPU time: " << time_sec << " s" << endl;
    cout << " Energy : " << m_energy << endl;
    cout << " Variance : " << m_variance << endl;
    cout << " Accepted step ratio : " << m_acceptRatio << endl;
    cout << endl;


}

void Sampler::computeAverages() {
    //Computing the averages of the sampled quantities.
    double steps_min_eq=m_system->getNumberOfMetropolisSteps()*(1-m_system->getEquilibrationFraction());
    m_energy = m_cumulativeEnergy / steps_min_eq;
    m_cumulativeEnergy2 =m_cumulativeEnergy2/ steps_min_eq;
    m_variance=m_cumulativeEnergy2-(m_energy*m_energy);
    m_acceptRatio = m_acceptedSteps / steps_min_eq;
    m_E_Lderiv=m_cumulativeE_Lderiv/steps_min_eq;
    m_E_Lderiv_expect=m_cumulativeE_Lderiv_expect/steps_min_eq;
}

double Sampler::computeVariance(std::vector<double> x_sample, double x_mean){
    double var_sum=0;
    for (int i=0; i<m_acceptedSteps; i++){
        var_sum+=pow((x_sample[i]-x_mean),2);
    }
    cout<<x_sample.size();
    return var_sum/(x_sample.size()-1);

}

//Here on its just functions for writing different quantities to file
//Benchmarking results
void Sampler::writeToFile(){
  ofstream myfile, myfiletime;
  string folderpart1, folderpart2;

  folderpart1 ="Results/"+m_system->getSampleMethod()+"/";

  int parti= m_system->getNumberOfParticles();
  int dimen= m_system->getNumberOfDimensions();

  std::string filename=folderpart1+"N="+std::to_string(parti)+"Dim="+std::to_string(dimen);
  std::string filenametime=folderpart1+"time/"+"N="+std::to_string(parti)+"Dim="+std::to_string(dimen);

  myfile.open(filename);
  myfiletime.open(filenametime);

  cout << "Mean energies are being written to file.."<<endl;
  myfiletime<<time_sec<<endl;
  myfiletime.close();

  for(int i=0; i<meanenergy_list.size(); i++){
    myfile<< std::fixed << std::setprecision(8) <<meanenergy_list[i]<<endl;
  }
  cout << "Done!"<<endl;
  cout<<endl;

  myfile.close();


}

/*
//Step sizes and time steps written to file
void Sampler::writeToFileSteps(std::vector<int> steps_list, std::vector<double> meanEL_list){
  ofstream myfile4;
  string folderpart1, folderpart2;
  
  if (m_system->getBruteforce()==true){
    folderpart1="Results/steps/bruteforce/";
    folderpart2="steplength"+std::to_string(m_system->getStepLength());
  }
  else {
    folderpart1 ="Results/steps/importancesampling/";
    folderpart2="timestep"+std::to_string(m_system->getTimeStep());
  }

  int parti= m_system->getNumberOfParticles();
  int dimen= m_system->getNumberOfDimensions();

  std::string filename=folderpart1+folderpart2+"N"+std::to_string(parti)+"Dim"+std::to_string(dimen)+".txt";
  myfile4.open(filename);
  cout << "Steps and energies are being written to file.."<<endl;
  for(int i=0; i<steps_list.size(); i++){
      myfile4<<steps_list[i]<<" "<<meanEL_list[i]<<endl;
  }
  cout << "Done!"<<endl;
  cout<<endl;

  myfile4.close();

}
*/