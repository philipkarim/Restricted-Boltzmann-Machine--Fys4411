#pragma once
#include <vector>
#include <Math/random.h>
#include <string>

using namespace std;

class System {
public:
    System();
    System(int seed);
    //Added aftertremoved with typos?
    /*
    void setPosition(const std::vector<double> &position);
    void adjustPosition(double change, int dimension);
    void setNumberOfDimensions(int numberOfDimensions);
    std::vector<double> getPosition() { return m_position; }
    */
    //
    bool metropolisStep             ();
    bool metropolisStepImportanceSampling();
    bool GibbsSampling();

    void runMetropolisSteps         (int RBMCycles, int numberOfMetropolisSteps);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    class Random*                   getRandomEngine()   { return m_random; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    int getRBMCycles()                  { return m_RBMCycles; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    
    void setSampleMethod                (char sampleMethod);
    char getSampleMethod()             {return m_sampleMethod;}

    void setTimeStep                    (double timeStep);
    double getTimeStep()                {return m_timeStep;}

    double getStepLength()               {return m_stepLength;}
    void setInteraction                 (bool interaction);
    double getInteraction()             {return m_interaction;}
    void setgeneralwtf                  (bool generalwtf);
    std::vector<double>get_energyarr()  { return m_energyarr; }

    //Try writing these in another way, maybe set dimension and particle from
    //main?
    void setNumberHiddenNodes       (int n_hidden){ m_numberHiddenNodes = n_hidden; }
    void setNumberVisibleNodes      (int n_visible){ m_numberVisibleNodes = n_visible; }
    void setNumberParticles(int numberOfParticles);
    void setNumberDimensions(int numberOfDimensions);



private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    int                             m_RBMCycles = 0;
    double                          m_equilibrationFraction = 0.0;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
    class Random*                   m_random = nullptr;
    std::vector<double>             m_energyarr;


   //Just some variables, mostly bools
    char m_sampleMethod;
    double m_stepLength=0.5;   //It said=0.1
    double m_timeStep=0.25;
    bool m_interaction;
    bool m_general_wtf;

    int                             m_numberHiddenNodes = 0;
    int                             m_numberVisibleNodes = 0;
};
