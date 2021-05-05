#include "randomuniform.h"
#include <iostream>
#include <cassert>
#include "../particle.h"
#include "../system.h"

using std::cout;
using std::endl;

RandomUniform::RandomUniform(System*    system,
                             int        numberOfDimensions,
                             int        numberOfParticles)  :
        InitialState(system) {
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;

    /* The Initial State class is in charge of everything to do with the
     * initialization of the system; this includes determining the number of
     * particles and the number of dimensions used. To make sure everything
     * works as intended, this information is passed to the system here.
     */
    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    setupInitialState();
}

void RandomUniform::setupInitialState() {
    m_x.zeros(m_nv);
    m_h.zeros(m_nh);
    m_a.zeros(m_nv);
    m_b.zeros(m_nh);
    m_w.zeros(m_nv, m_nh);

    std::uniform_real_distribution<double> uniform_weights(-m_initialization, m_initialization);
    std::uniform_real_distribution<double> uniform_position(-0.5, 0.5);
    std::normal_distribution<double> normal_weights(0, m_initialization);

    // initialize weights according to either a uniform or gaussian distribution
    if (m_normaldistr == 0){
        for (int i=0; i<m_nv; i++){
            m_a[i] = uniform_weights(m_randomEngine);
            for (int j=0; j<m_nh; j++){
                m_w(i, j) = uniform_weights(m_randomEngine);
            }
        }

        for (int j=0; j<m_nh; j++){
            m_b[j] = uniform_weights(m_randomEngine);
        }
    }

    else if (m_normaldistr == 1){
        for (int i=0; i<m_nv; i++){
            m_a[i] = normal_weights(m_randomEngine);
            for (int j=0; j<m_nh; j++){
                m_w(i, j) = normal_weights(m_randomEngine);
            }
        }

        for (int j=0; j<m_nh; j++){
            m_b[j] = normal_weights(m_randomEngine);
        }
    }

    for (int i=0; i<m_nv; i++){
        m_x[i] = uniform_position(m_randomEngine);
    }

}