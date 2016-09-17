#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace std;

int main()
{
    int nParticles = 5;
    int SPB_dim = 10;                       //dimension (dim) of single-particle basis (SPB)
    Eigen::MatrixXd C(SPB_dim, SPB_dim);    //HF coefficients

    Eigen::MatrixXf densityMat;
    densityMat.setIdentity(SPB_dim, SPB_dim);

    double sum = 0;
    for (int gamma=0; gamma<SPB_dim; gamma++){
        for (int delta=0; delta<SPB_dim; delta++){
            sum = 0;
            for (int i=0; i<nParticles; i++){
                sum += C(gamma,i)*C(delta,i);
            }
            densityMat(gamma,delta) = sum;
        }
    }

    int maxHFiter       = 100;
    double eps          = 1e-5;
    double difference   = 1;

    int hf_count = 0;
    Eigen::Matrix3d oldEnergies;
    Eigen::Matrix3d newEnergies;

    while (hf_count < maxHFiter && difference > eps){
        Eigen::MatrixXd HFmatrix(SPB_dim, SPB_dim);

        for (int alpha=0; alpha<SPB_dim; alpha++){
            for (int beta=0; beta<SPB_dim; beta++){
                double sumFockTerm = 0;
                for (int gamma=0; gamma<SPB_dim; gamma++){
                    for (int delta=0; delta<SPB_dim; delta++){
                        //sumFockTerm += densityMat(gamma,delta)*twoBodyInteraction(alpha,gamma,beta,delta);
                    }
                }
                HFmatrix(alpha,beta) = sumFockTerm;
                if (beta == alpha){
                   // HFmatrix(alpha,beta) = oneBodyInteraction(alpha); //one-body Hamiltonian
                }
            }
        }


        Eigen::EigenSolver<Eigen::MatrixXd> es(HFmatrix);     //creates object "es" consisting of eigenvals/vecs of HFmatrix
        Eigen::VectorXd spEnergies = es.eigenvalues();
        C = es.eigenvectors();

        Eigen::MatrixXf densityMat(SPB_dim, SPB_dim);

        double sum = 0;
        for (int gamma=0; gamma<SPB_dim; gamma++){
            for (int delta=0; delta<SPB_dim; delta++){
                sum = 0;
                for (int i=0; i<nParticles; i++){
                    sum += C(gamma,i)*C(delta,i);
                }
                densityMat(gamma,delta) = sum;
            }
        }

        newEnergies = spEnergies;

        sum = 0;
        for (int i=0; i<SPB_dim; i++){
            sum += abs(newEnergies(i) - oldEnergies(i))/SPB_dim;
        }
        difference = sum;
        oldEnergies = newEnergies;

        cout << "Single-particle energies, ordering may have changed " << endl;
        for (int i=0; i<SPB_dim; i++){
            cout << "i,     oldEnergies:", oldEnergies(i) << endl;
        }
        hf_count += 1;
    }
}

double oneBodyInteraction(alpha){
    return spEnergies(alpha);
}

double twoBodyInteraction(alpha,gamma,beta,delta){

}

//setUpBasis takes in string "type", spesifying what inital basis to use.
//For nuclear matter, use plane wave basis
double setUpBasis(type,){

}
