//imported files
//other libraries
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <chrono>
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

int main()
{
    double  eps     =   1e-15;                      //remember to adjust setprecision in master when changing this
    double  conFac  =   1;                          //convergence factor

    bool    intermediates = true;                   //turn on/off intermediates in CCD eqs
    bool    CCDT          = true;                  //turn on/off CCDT-1
    bool    timer         = true;                   //turn on/off timer

    bool    makeData      = false;                  //choose to write to file for a range of shells

    //Eigen::initParallel();

    if (makeData == false){
        //we use natural units
        const double  pi      =   M_PI;
        int     Nh      =   14;							//number of particles
        int     Nb      =   3;							//number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 can at max have N=14)
        double  rs      =   1;                          //Wigner Seitz radius
        double  rb      =   1;                          //Bohr radius [MeV^-1]
        double  m       =   1;//                //electron mass [MeV] (1 for HEG, 939.565 for MP)
        //double  m       =   939.565;
        //double  rho     =   0.2;
        double  L3      =   4*pi*Nh*rs/3;               //box volume
        //double  L3      =   Nh/rho;
        double  L2      =   pow(L3, 2./3.);
        double  L1      =   pow(L3, 1./3.);

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

    }
    else if (makeData == true){

        //we use natural units
        const double  pi      =   M_PI;

        int lower_bound = 2; int upper_bound = 41;  //set lower and upper limits on shells to be used

        int     Nh      =   14;							//number of particles
        double  rs      =   1;                          //Wigner Seitz radius
        double  rb      =   1;                          //Bohr radius [MeV^-1]
        //double  m       =   1;//                //electron mass [MeV] (1 for HEG, 939.565 for MP)
        double  m       =   939.565;
        double  rho     =   0.2;
        //double  L3      =   4*pi*Nh*rs/3;               //box volume
        double  L3      =   Nh/rho;
        double  L2      =   pow(L3, 2./3.);
        double  L1      =   pow(L3, 1./3.);

        double ECC;
        ofstream myfile;

        myfile.open("Nh_14_MP_rho02.txt");
        for (int Nb=lower_bound; Nb<upper_bound; Nb++){
            Master* master = new Master;
            master->setSize(Nh, Nb);

            master->setSystem(new MP(master, m, L3, L2, L1));

            master->setTriples(CCDT);
            master->setIntermediates(intermediates);
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

            myfile <<  Nb << " " << "&" << " " << master->m_Ns << " " << "&" << " " << std::setprecision(16) << ECC << " \\\\ " << "\n";

        }
        myfile.close();
    }

    return 0;

}
