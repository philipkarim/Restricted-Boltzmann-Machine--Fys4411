#pragma once
#include "initialstate.h"
#include <armadillo>

class RandomUniform : public InitialState {
public:
    RandomUniform(System* system, int n_hidden, 
                    int n_visible, bool gaussian, 
                    double initialization);
    void setupInitialState();

private:    
    double m_nv;
    double m_nh;
    double m_initialization;
    double m_sigma;
    bool m_normaldistr;


};

