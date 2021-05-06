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

    arma::vec initial_x;      // visible nodes (i.e. position)
    arma::vec initial_h;      // hidden nodes
    arma::vec initial_a;      // visible bias
    arma::vec initial_b;      // hidden bias
    //This is actually a matrix
    arma::mat initial_w;      // interaction of biases

};

