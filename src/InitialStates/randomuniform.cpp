#include "randomuniform.h"
#include "initialstate.h"
#include <iostream>
#include <cassert>
#include "../particle.h"
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

//#include <Eigen/Dense>
#include<armadillo>
using namespace std;
using namespace arma;
//using namespace Eigen
RandomUniform::RandomUniform(System* system, int n_hidden, 
                            int n_visible, bool gaussian, 
                            double initialization):
        InitialState(system) {
    assert(n_hidden >= 0 && n_visible >= 0);
    m_nh = n_hidden;
    m_nv = n_visible;
    m_normaldistr = gaussian;
    m_initialization = initialization;

    /* The Initial State class is in charge of everything to do with the
     * initialization of the system; this includes determining the number of
     * nodes and the number of the distribution used. To make sure everything
     * works as intended, this information is passed to the system here.
     */
    m_system->setNumberOfHN(n_hidden);
    m_system->setNumberOfVN(n_visible);

    setupInitialState();
}

void RandomUniform::setupInitialState() {
    vec initial_x(m_nv);      // visible nodes (i.e. position)
    vec initial_h(m_nh);      // hidden nodes
    vec initial_a(m_nv);      // visible bias
    vec initial_b(m_nh);      // hidden bias
    //This is actually a matrix
    mat initial_w(m_nv, m_nh);      // interaction of biases

    initial_x.zeros();
    initial_h.zeros();
    initial_a.zeros();
    initial_b.zeros();
    initial_w.zeros();

    random_device rd;
    mt19937_64 gen(rd());

    // Set up the distribution for x \in [[x, x],(can use multiple configurations)
    uniform_real_distribution<double> UniformNumberGenerator(-0.5,0.5);
    uniform_real_distribution<double> uniform_weights(-m_initialization, m_initialization);
    normal_distribution<double> normal_weights(0, m_initialization);

    // initialize weights according to either a uniform or gaussian distribution
    if (m_normaldistr == 0){
        for (int i=0; i<m_nv; i++){
            initial_a[i] = uniform_weights(gen);
            for (int j=0; j<m_nh; j++){
                initial_w(i, j) = uniform_weights(gen);
            }
        }

        for (int j=0; j<m_nh; j++){
            initial_b[j] = uniform_weights(gen);
        }
    }

    else if (m_normaldistr == 1){
        for (int i=0; i<m_nv; i++){
            initial_a[i] = normal_weights(gen);
            for (int j=0; j<m_nh; j++){
                initial_w(i, j) = normal_weights(gen);
            }
        }

        for (int j=0; j<m_nh; j++){
            initial_b[j] = normal_weights(gen);
        }
    }

    for (int i=0; i<m_nv; i++){
        initial_x[i] = UniformNumberGenerator(gen);
    }

    m_system->getWaveFunction()->set_a(initial_a);
    m_system->getWaveFunction()->set_b(initial_b);
    m_system->getWaveFunction()->set_w(initial_w);
    m_system->getWaveFunction()->set_X(initial_x);
    m_system->getWaveFunction()->set_h(initial_h);

}