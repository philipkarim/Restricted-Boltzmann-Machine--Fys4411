#pragma once
#include <vector>

using namespace std;
class InitialState {
public:
    InitialState(class System* system);
    virtual void setupInitialState() = 0;

protected:
    class System* m_system = nullptr;
    vector<double> initial_x;      // visible nodes (i.e. position)
    vector<double> initial_h;      // hidden nodes
    vector<double> initial_a;      // visible bias
    vector<double> initial_b;      // hidden bias
    //This is actually a matrix
    vector<vector<double>> initial_w;      // interaction of biases

};

