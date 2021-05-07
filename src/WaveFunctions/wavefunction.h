#pragma once
#include <iostream>
#include <armadillo>



class WaveFunction {
public:
    WaveFunction(class System* system);

    virtual double evaluate(arma::vec position)=0;
    //void setupInitialState();
    virtual double computeDoubleDerivative(arma::vec position)=0;
    virtual double computeDerivative(arma::vec position)=0;
    virtual double sigmoid(double x)=0;
    virtual double v(int j)=0;


    virtual void set_X(arma::vec X)=0;
    virtual void set_h(arma::vec h)=0;
    virtual void set_a(arma::vec a)=0;
    virtual void set_b(arma::vec b)=0;
    virtual void set_w(arma::mat w)=0;

    virtual arma::vec get_X()=0;
    virtual arma::vec get_h()=0;
    virtual arma::vec get_a()=0;
    virtual arma::vec get_b()=0;
    //This is actually a matrix. Use armadillo or eigen maybe?
    virtual arma::mat get_w()=0;

    //virtual void setupInitialState();


//private:
    //Not sure about this one, but move it in the other too up or down a class
  //  void setupInitialState() = 0;

protected:
    arma::vec m_x;      // visible nodes (i.e. position)
    arma::vec m_h;      // hidden nodes
    arma::vec m_a;      // visible bias
    arma::vec m_b;      // hidden bias
    //This is actually a matrix
    arma::mat m_w;      // interaction of biases

    double m_nv=0;
    double m_nh=0;
    //Maybe remove these two_?
    int m_part = 0;     // number of particles
    int m_dim = 0;      // number of dimensions

    class System* m_system = nullptr;

};
