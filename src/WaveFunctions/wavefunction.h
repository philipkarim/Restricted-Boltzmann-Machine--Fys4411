#pragma once
#include <iostream>
#include <armadillo>

using namespace arma;

class WaveFunction {
public:
    WaveFunction(class System* system);

    virtual double evaluate(vec position)=0;
    //void setupInitialState();
    virtual double computeDoubleDerivative(vec position)=0;
    virtual double computeDerivative(vec position)=0;
    //virtual double sigmoid(double x)=0;
    //virtual double sigmoid_input(int x)=0;
    virtual double computeQuantumForce(vec position, int index)=0;


    virtual void set_X(vec X)=0;
    virtual void set_h(vec h)=0;
    virtual void set_a(vec a)=0;
    virtual void set_b(vec b)=0;
    virtual void set_w(mat w)=0;

    virtual vec get_X()=0;
    virtual vec get_h()=0;
    virtual vec get_a()=0;
    virtual vec get_b()=0;
    //This is actually a matrix. Use armadillo or eigen maybe?
    virtual mat get_w()=0;
    virtual double getSigma()=0;

    //virtual void setupInitialState();


//private:
    //Not sure about this one, but move it in the other too up or down a class
  //  void setupInitialState() = 0;

protected:
    vec m_x;      // visible nodes (i.e. position)
    vec m_h;      // hidden nodes
    vec m_a;      // visible bias
    vec m_b;      // hidden bias
    //This is actually a matrix
    mat m_w;      // interaction of biases

    double m_nv=0;
    double m_nh=0;
    //Maybe remove these two_?
    int m_part = 0;     // number of particles
    int m_dim = 0;      // number of dimensions

    class System* m_system = nullptr;
private:
    virtual double sigmoid(double x)=0;
    virtual double sigmoid_input(int x)=0;

};
