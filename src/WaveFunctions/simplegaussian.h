#pragma once
#include "wavefunction.h"
#include <iostream>
#include <vector>

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, int n_hidden, int n_visible, int part, int dim, double sigma, bool gaussian, double initialization);
    double evaluate(std::vector<double> position);
    //double setupInitialState();
    double computeDoubleDerivative(std::vector<double> position);
    double computeDerivative(std::vector<double> position);

    std::vector<double> set_X(std::vector<double> X){ m_x = X; }
    std::vector<double> set_h(std::vector<double> h){ m_h = h; }
    std::vector<double> set_a(std::vector<double> a){ m_a = a; }
    std::vector<double> set_b(std::vector<double> b){ m_b = b; }
    //This is also a matrix 
    std::vector<double> set_w(std::vector<double> w) {m_w = w; }

    std::vector<double> get_X(){ return m_x; }
    std::vector<double> get_h(){ return m_h; }
    std::vector<double> get_a(){ return m_a; }
    std::vector<double> get_b(){ return m_b; }
    //This is actually a matrix. Use armadillo or eigen maybe?
    std::vector<double> get_w(){ return m_w; }
    
    void setupInitialState();


private:
    std::vector<double> m_x;      // visible nodes (i.e. position)
    std::vector<double> m_h;      // hidden nodes
    std::vector<double> m_a;      // visible bias
    std::vector<double> m_b;      // hidden bias
    //This is actually a matrix
    std::vector<double> m_w;      // interaction of biases

    double m_nv;
    double m_nh;
    double m_initialization;
    double m_sigma;
    string m_gaussian;

    //Should this be here?
    //void setupInitialState();
};
