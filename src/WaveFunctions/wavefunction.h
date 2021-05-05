#pragma once
#include <iostream>
#include <vector>

using namespace std;

class WaveFunction {
public:
    WaveFunction(class System* system);

    double evaluate(vector<double> position);
    //void setupInitialState();
    double computeDoubleDerivative(vector<double> position);
    double computeDerivative(vector<double> position);

    virtual void set_X(vector<double> X)=0;
    virtual void set_h(vector<double> h)=0;
    virtual void set_a(vector<double> a)=0;
    virtual void set_b(vector<double> b)=0;
    virtual void set_w(vector<vector<double>> w)=0;

    virtual vector<double> get_X()=0;
    virtual vector<double> get_h()=0;
    virtual vector<double> get_a()=0;
    virtual vector<double> get_b()=0;
    //This is actually a matrix. Use armadillo or eigen maybe?
    virtual vector<vector<double>> get_w()=0;

    //virtual void setupInitialState();


//private:
    //Not sure about this one, but move it in the other too up or down a class
  //  void setupInitialState() = 0;

protected:
    vector<double> m_x;      // visible nodes (i.e. position)
    vector<double> m_h;      // hidden nodes
    vector<double> m_a;      // visible bias
    vector<double> m_b;      // hidden bias
    //This is actually a matrix
    vector<vector<double>> m_w;      // interaction of biases

    double m_nv=0;
    double m_nh=0;
    //Maybe remove these two_?
    int m_part = 0;     // number of particles
    int m_dim = 0;      // number of dimensions

    class System* m_system = nullptr;

};
