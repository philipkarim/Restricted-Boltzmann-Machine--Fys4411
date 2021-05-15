#include "neuralstate.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>
#include <armadillo>

using namespace arma;
//using namespace std;
NeuralState::NeuralState(System* system, int part, 
                                int dim, double sigma) :
        WaveFunction(system) {
    m_sigma = sigma;
    m_system->setNumberParticles(part);
    m_system->setNumberDimensions(dim);

}

double NeuralState::evaluate(vec position) {
    //Implementation of wavefunction at the given position
    double exponent_one=0;
    double product_term=1;
    double sum_in_product=0;
    double psi_value;

    for (int i=0; i<m_system->getNumberOfVN(); i++){
        exponent_one+=((position[i]-m_a[i])*(position[i]-m_a[i]))/(2*m_sigma*m_sigma);
    }

    //Might have to transpose some of the matrixes and vectors
    for (int j=0; j<m_system->getNumberOfHN(); j++){
        for (int ii=0; ii<m_system->getNumberOfVN(); ii++){
            sum_in_product+=(position[ii]*m_w[ii, j])/(m_sigma*m_sigma);
            product_term*=(1+exp(m_b[j]+sum_in_product));
        }
    }
    psi_value=exp(-exponent_one)*product_term;

    return psi_value;
}



    //Return a double value
    //return 1.};

double NeuralState::computeDoubleDerivative(vec position) {
    //Computes the value of the analytical double derivative for the non interacting case.
    double sum_M=0;
    double sum_N=0;
    double sig_inp;
    for (int i=0; i<m_system->getNumberOfVN(); i++){
        sum_M-=1/(m_sigma*m_sigma);
        for (int j=0; j<m_system->getNumberOfHN(); j++){
            sig_inp=sigmoid_input(j);
            sum_N+=(m_w(i,j)*m_w(i,j))/(pow(m_sigma, 4))*sigmoid(sig_inp)*sigmoid(-sig_inp);
        }
        sum_M+=sum_N;
    }

    //cout<<"________________\n"<<(get_w())<<"____________";
    return sum_M;
}

/*
    int M = m_system->getNumberOfVN();
    int N = m_system->getNumberOfHN();

    arma::vec S; S.zeros(N);
    arma::vec S_tilde; S_tilde.zeros(N);
    arma::vec one_vector; one_vector.ones(M);

    for (int j = 0; j<N; j++){
        S[j] = sigmoid(v(j));
        S_tilde[j] = sigmoid(v(j))*sigmoid(-v(j));
    }

    arma::vec O = m_b + (m_x.t()*m_w).t()/(m_sigma*m_sigma);

    double Ek = 0.0;
    double term1 = 0.0;
    double term2 = 0.0;
    double term3 = 0.0;
    double term4 = 0.0;
    for (int i=0; i<m_nv; i++){
        double sum1 = 0.0;
        double sum2 = 0.0;
        double sum3 = 0.0;
        for (int j=0; j<m_nh; j++){
            sum1 += m_w(i, j)*S[j];
            sum2 += m_w(i, j)*S[j];
            sum3 += m_w(i, j)*m_w(i, j)*S_tilde[j];
        }
        term1 += (m_x[i] - m_a[i])*(m_x[i] - m_a[i]);
        term2 += -2*(m_x[i] - m_a[i])*sum1;
        term3 += sum2*sum2;
        term4 += sum3;
    }

    Ek = -1.0/(2*pow(m_sigma, 4))*(term1 + term2 + term3 + term4) + m_nv/(2*m_sigma*m_sigma);

    cout<<Ek;
    return Ek;
  }
  */

//Computes the value of the analytically first derivative 
double NeuralState::computeDerivative(vec position) {
    double first_sum=0;
    double sec_sum=0;
    
    for (int i =0; i<m_system->getNumberOfVN(); i++){
        first_sum-=(position[i]-m_a[i])/(m_sigma*m_sigma);
        for (int j=0; j<m_system->getNumberOfHN(); j++){
            sec_sum+=m_w(i,j)/(m_sigma*m_sigma)*sigmoid(sigmoid_input(j));
        }
        first_sum+=sec_sum;
    }
    //Return a double value
    return first_sum;
  }


//Just the sigmoid function to be used in the derivative functions
double NeuralState::sigmoid(double x){
    return (1/(1+exp(-x)));
}

//Computes the input of the sigmoid function
double NeuralState::sigmoid_input(int x){
    double sum=1.;
    for (int i=0; i<m_system->getNumberOfVN(); i++){
        sum+=m_x(i)*m_w(i,x)/(m_sigma*m_sigma);
    }
    //return 1;
    return m_b(x)+sum;
}

double NeuralState::computeQuantumForce(vec position, int index){

    arma::vec O = m_b + (m_x.t()*m_w).t()/(m_sigma*m_sigma);

    double deriv=0.0;
    double sum1 = 0.0;
    for(int j=0; j<m_nh; j++){
        sum1 += m_w(index, j)/(m_sigma*m_sigma*(1.0+exp(-O[j])));
    }
    deriv = -(position[index]-m_a[index])/(m_sigma*m_sigma) + sum1/(m_sigma*m_sigma);

    return 2*deriv;
}