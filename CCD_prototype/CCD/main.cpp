//imported files
//other libraries
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <chrono>
#include <mpi.h>
#include <omp.h>
#include <eigen3/Eigen/Dense>
//#include <eigen3/Eigen/Core>

typedef std::chrono::high_resolution_clock Clock;   //needed for timing

//author made files
#include "makestatespace.h"
#include "master.h"
#include "Systems/system.h"
#include "Systems/heg.h"
#include "Systems/mp.h"
#include "makeampmat.h"
#include "makeintmat.h"


using namespace std;

int main(int argc, char** argv)
{

    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); //sets world_size equal to world size
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); //sets world_rank equal to world rank

    MPI_Status status;

    double  eps     =   1e-16;                      //remember to adjust setprecision in master when changing this
    double  conFac  =   1;                          //convergence factor

    bool    intermediates = true;                   //turn on/off intermediates in CCD eqs
    bool    CCDT          = true;                   //turn on/off CCDT-1
    bool    timer         = true;                   //turn on/off timer
    bool    relaxation    = true;                  //turn on/off relaxation when updating amplitudes
    double  alpha         = 0.8;                    //relaxation parameter

    bool    makeData      = false;                  //choose to write to file for a range of shells

    Eigen::initParallel();

    int threads = 2;
    omp_set_num_threads(threads);
    Eigen::setNbThreads(threads);

    int n = Eigen::nbThreads( );

    std::cout << n << std::endl;

    if (makeData == false){
        //we use natural units
        const double  pi      =   M_PI;
        int     Nh      =   14;							//number of particles
        int     Nb      =   4;							//number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 can at max have N=14)
        double  rs      =   1;                          //Wigner Seitz radius
        double  rb      =   1.;                          //Bohr radius [MeV^-1]
        double  m       =   1.;//                //electron mass [MeV] (1 for HEG, 939.565 for MP)
        //double  m       =   939.565;
        //double  rho     =   0.2;
        double  r1      =   pow(rs*rb, 3);
        double  L3      =   4.*pi*Nh*r1/3.;               //box volume
        //double  L3      =   Nh/rho;
        double  L2      =   pow(L3, 2./3.);
        double  L1      =   pow(L3, 1./3.);
        //cout << L1 << endl;

        Master* master = new Master;
        master->setSize(Nh, Nb);

        master->setSystem(new HEG(master, m, L3, L2, L1));

        master->setTriples(CCDT);
        master->setIntermediates(intermediates);
        master->setRelaxation(relaxation, alpha);
        master->setTimer(timer);


        cout << "C++ code" << endl;

        auto t1 = Clock::now();

        if (CCDT){
            double ECCDT = master->CC_master(eps, conFac);
            cout << "Delta ECCDT-1: "<< ECCDT << endl;
        }
        else{
            double ECCD = master->CC_master(eps, conFac);
            cout << "Delta ECCD: "<< ECCD << endl;
        }

        auto t2 = Clock::now();

        if (intermediates){
            std::cout << "Total time used: "
                      << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
                      << " seconds, with intermediates ON" << std::endl;
        }
        else{
            std::cout << "Total time used: "
                      << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
                      << " seconds, with intermediates OFF" << std::endl;
        }

    }
    else if (makeData == true){

        //we use natural units
        const double  pi      =   M_PI;

        int lower_bound = 2; int upper_bound = 34;  //set lower and upper limits on shells to be used

        //int     Nb      =   3;
        int     Nh      =   14;							//number of particles
        //double  rs      =   1;                          //Wigner Seitz radius
        //double  rb      =   1;                          //Bohr radius [MeV^-1]
        //double  m       =   1;//                //electron mass [MeV] (1 for HEG, 939.565 for MP)
        //double  m       =   939.565;
        //double  rho     =   0.5;
        //double  L3      =   4*pi*Nh*rs/3;               //box volume
        //double  L3      =   Nh/rho;
        //double  L2      =   pow(L3, 2./3.);
        //double  L1      =   pow(L3, 1./3.);


        int prev_num_states = 0;
        int curr_num_states;

        for (int j=2; j<20; j++){
            double ECC;
            ofstream myfile;
            ostringstream os;
        double  rs      =   j/(double)10;                          //Wigner Seitz radius
        os << "Nh_14_HEG_rs" << rs*10 <<".txt";
        string s = os.str();
        myfile.open(s);
        for (int i=lower_bound; i<upper_bound; i++){
            int Nb = i;
            double  rb      =   1;
            //std::cout << rs << std::endl;
            double  m       =   1;                //electron mass [MeV] (1 for HEG, 939.565 for MP)
            double  r1      =   pow(rs*rb, 3);
            double  L3      =   4.*pi*Nh*r1/3.;               //box volume
            double  L2      =   pow(L3, 2./3.);
            double  L1      =   pow(L3, 1./3.);

            Master* master = new Master;
            master->setSize(Nh, Nb);

            master->setSystem(new HEG(master, m, L3, L2, L1));

            curr_num_states = master->m_Ns;
            //std::cout << curr_num_states << std::endl;

            if (curr_num_states != prev_num_states){

            prev_num_states = curr_num_states;

            master->setTriples(CCDT);
            master->setIntermediates(intermediates);
            master->setRelaxation(relaxation, alpha);
            master->setTimer(timer);


            cout << "C++ code" << endl;

            auto t1 = Clock::now();

            if (CCDT){
                ECC = master->CC_master(eps, conFac);
                cout << "Delta ECCDT-1: "<< ECC << endl;
            }
            else{
                ECC = master->CC_master(eps, conFac);
                cout << "Delta ECCD: "<< ECC << endl;
            }

            auto t2 = Clock::now();

            if (intermediates){
                std::cout << "Total time used: "
                          << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
                          << " seconds, with intermediates ON" << std::endl;
            }
            else{
                std::cout << "Total time used: "
                          << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
                          << " seconds, with intermediates OFF" << std::endl;
            }

            myfile <<  Nb << " " << "&" << " " << master->m_Ns << " " << "&" << " " << std::setprecision(16) << ECC << " \\\\ " << "\n" << std::flush;
            }

        }
        myfile.close();
    }
    }

    MPI_Finalize();

    return 0;

}
