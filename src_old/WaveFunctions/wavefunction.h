#pragma once
#include <iostream>
#include <vector>

class WaveFunction {
public:
    WaveFunction(class System* system);

    double evaluate(std::vector<double> position)=0;
    double setupInitialState()=0;
    double computeDoubleDerivative(std::vector<double> position)=0;
    double computeDerivative(std::vector<double> position)=0;

    virtual std::vector<double> set_X(std::vector<double> X)=0;
    virtual std::vector<double> set_h(std::vector<double> h)=0;
    virtual std::vector<double> set_a(std::vector<double> a)=0;
    virtual std::vector<double> set_b(std::vector<double> b)=0;
    virtual std::vector<double> set_w(std::vector<double> w)=0;

    virtual std::vector<double> get_X()=0;
    virtual std::vector<double> get_h()=0;
    virtual std::vector<double> get_a()=0;
    virtual std::vector<double> get_b()=0;
    //This is actually a matrix. Use armadillo or eigen maybe?
    virtual std::vector<double> get_w()=0;

    virtual void setupInitialState() = 0;


//private:
    //Not sure about this one, but move it in the other too up or down a class
  //  void setupInitialState() = 0;

protected:
    std::vector<double> m_x;      // visible nodes (i.e. position)
    std::vector<double> m_h;      // hidden nodes
    std::vector<double> m_a;      // visible bias
    std::vector<double> m_b;      // hidden bias
    //This is actually a matrix
    std::vector<double> m_w;      // interaction of biases

    double m_nv=0;
    double m_nh=0;
    //Maybe remove these two_?
    int m_part = 0;     // number of particles
    int m_dim = 0;      // number of dimensions

    class System* m_system = nullptr;

};
