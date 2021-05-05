#pragma once
#include "wavefunction.h"
#include <iostream>
#include <vector>

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, int part, int dim, double sigma);
    double evaluate(std::vector<double> position);
    double computeDoubleDerivative(std::vector<double> position);
    double computeDerivative(std::vector<double> position);

    void set_X(std::vector<double> X){ m_x = X; }
    void set_h(std::vector<double> h){ m_h = h; }
    void set_a(std::vector<double> a){ m_a = a; }
    void set_b(std::vector<double> b){ m_b = b; }
    //This is also a matrix 
    void set_w(std::vector<vector<double>> w) {m_w = w; }

    std::vector<double> get_X(){ return m_x; }
    std::vector<double> get_h(){ return m_h; }
    std::vector<double> get_a(){ return m_a; }
    std::vector<double> get_b(){ return m_b; }
    //This is actually a matrix. Use armadillo or eigen maybe?
    std::vector<vector<double>>get_w(){ return m_w; }
    
    //void setupInitialState();

/*
private:

    std::vector<double> m_x;      // visible nodes (i.e. position)
    std::vector<double> m_h;      // hidden nodes
    std::vector<double> m_a;      // visible bias
    std::vector<double> m_b;      // hidden bias
    //This is actually a matrix
    std::vector<vector<double>> m_w;      // interaction of biases
*/
private:
    double m_sigma;
    std::vector<double> m_x;      // visible nodes (i.e. position)
    std::vector<double> m_h;      // hidden nodes
    std::vector<double> m_a;      // visible bias
    std::vector<double> m_b;      // hidden bias
    //This is actually a matrix
    std::vector<vector<double>> m_w;      // interaction of biases
    /*
    double m_nv;
    double m_nh;
    double m_initialization;
    
    */
    //Should this be here?
    //void setupInitialState();
};
