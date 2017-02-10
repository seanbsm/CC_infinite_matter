//imported files
//other libraries
#include <iostream>
#include <math.h>
#include <chrono>
#include <eigen3/Eigen/Dense>

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

int main()
{
//we use natural units
    const double  pi      =   M_PI;
    int     Nh      =   14;							//number of particles
    int     Nb      =   3;							//number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 can at max have N=14)
    double  rs      =   1;                          //Wigner Seitz radius
    double  rb      =   1;                          //Bohr radius [MeV^-1]
    double  m       =   1;//939.565;                    //electron mass [MeV] (1 for HEG, 939.565 for MP)
    //double  rho     =   0.0001;
    double  L3      =   4*pi*Nh*rs/3;               //box volume
    //double  L3      =   Nh/rho;
    double  L2      =   pow(L3, 2./3);
    double  L1      =   pow(L3, 1./3);

    double  eps     =   1e-12;
    double  conFac  =   1;                          //convergence factor

    bool    intermediates = true;                   //turn on/off intermediates in CCD eqs
    bool    CCDT          = true;                  //turn on/off CCDT-1
    bool    timer         = true;                   //turn on/off timer



    Master* master = new Master;
    master->setSize(Nh, Nb);

    master->setSystem(new HEG(master, m, L3, L2, L1));

    master->setTriples(CCDT);
    master->setIntermediates(intermediates);
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




    return 0;

    /*
    int lower_bound = 3;
    int upper_bound = 30;
    std::vector<double> Es;
    for (int Nb=lower_bound; Nb<upper_bound; Nb++){
        Master* master = new Master;
        master->setSize(Nh, Nb);

        master->setSystem(new HEG(master, m, L3, L2, L1));

        master->setTriples(CCDT);
        master->setIntermediates(intermediates);
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
            Es.push_back(ECCD);
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

    for (int i=0; i<Es.size(); i++){
        cout << Es[i] << endl;
    }

    return 0;
    */
}
