#pragma once
#include "wavefunction.h"
#include <iostream>
#include <armadillo>

class NeuralState : public WaveFunction {
public:
    NeuralState(class System* system, int part, int dim, double sigma);
    double evaluate(arma::vec position);
    double computeDoubleDerivative(arma::vec position);
    double computeDerivative(arma::vec position);

    double sigmoid(double x);
    double sigmoid_input(int x);

    void set_X(arma::vec X){ m_x = X; }
    void set_h(arma::vec h){ m_h = h; }
    void set_a(arma::vec a){ m_a = a; }
    void set_b(arma::vec b){ m_b = b; }
    //This is also a matrix 
    void set_w(arma::mat w) {m_w = w; }

    arma::vec get_X(){ return m_x; }
    arma::vec get_h(){ return m_h; }
    arma::vec get_a(){ return m_a; }
    arma::vec get_b(){ return m_b; }
    //This is actually a matrix. Use armadillo or eigen maybe?
    arma::mat get_w(){ return m_w; }
    
    //void setupInitialState();

/*
private:

    arma::vec m_x;      // visible nodes (i.e. position)
    arma::vec m_h;      // hidden nodes
    arma::vec m_a;      // visible bias
    arma::vec m_b;      // hidden bias
    //This is actually a matrix
    std::vector<vector<double>> m_w;      // interaction of biases
*/
private:
    double m_sigma;
    arma::vec m_x;      // visible nodes (i.e. position)
    arma::vec m_h;      // hidden nodes
    arma::vec m_a;      // visible bias
    arma::vec m_b;      // hidden bias
    //This is actually a matrix
    arma::mat m_w;      // interaction of biases
    /*
    double m_nv;
    double m_nh;
    double m_initialization;
    
    */
    //Should this be here?
    //void setupInitialState();
};
