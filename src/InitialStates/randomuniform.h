#pragma once
#include "initialstate.h"
#include <vector>

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

    std::vector<double> initial_x;      // visible nodes (i.e. position)
    std::vector<double> initial_h;      // hidden nodes
    std::vector<double> initial_a;      // visible bias
    std::vector<double> initial_b;      // hidden bias
    //This is actually a matrix
    std::vector<vector<double>> initial_w;      // interaction of biases

};

