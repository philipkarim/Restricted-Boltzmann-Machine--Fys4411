#pragma once
#include "wavefunction.h"
#include <iostream>
#include <armadillo>

using namespace arma;

class NeuralState : public WaveFunction {
public:
    NeuralState(class System* system, int part, int dim, double sigma);
    double evaluate(vec position);
    double computeDoubleDerivative(vec position);
    double computeDerivative(vec position);

    //double sigmoid(double x);
    //double sigmoid_input(int x);

    void set_X(vec X){ m_x = X; }
    void set_h(vec h){ m_h = h; }
    void set_a(vec a){ m_a = a; }
    void set_b(vec b){ m_b = b; }
    //This is also a matrix 
    void set_w(mat w) {m_w = w; }

    vec get_X(){ return m_x; }
    vec get_h(){ return m_h; }
    vec get_a(){ return m_a; }
    vec get_b(){ return m_b; }
    //This is actually a matrix. Use armadillo or eigen maybe?
    mat get_w(){ return m_w; }
    
    double getSigma(){return m_sigma;}
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
    vec m_x;      // visible nodes (i.e. position)
    vec m_h;      // hidden nodes
    vec m_a;      // visible bias
    vec m_b;      // hidden bias
    //This is actually a matrix
    mat m_w;      // interaction of biases
    /*
    double m_nv;
    double m_nh;
    double m_initialization;
    
    */
    //Should this be here?
    //void setupInitialState();
    double sigmoid(double x);
    double sigmoid_input(int x);
};
