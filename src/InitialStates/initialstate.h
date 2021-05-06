#pragma once
#include <armadillo>

using namespace std;
class InitialState {
public:
    InitialState(class System* system);
    virtual void setupInitialState() = 0;

protected:
    class System* m_system = nullptr;
    //virtual arma::vec initial_x=0;      // visible nodes (i.e. position)
    //virtual arma::vec initial_h;      // hidden nodes
    //virtual arma::vec initial_a;      // visible bias
    //virtual arma::vec initial_b;      // hidden bias
    //This is actually a matrix
    //virtual arma::mat initial_w;      // interaction of biases

};

