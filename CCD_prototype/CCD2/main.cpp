//imported files
//other libraries
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <fstream>
#include <math.h>
#include <chrono>
#include <mpi.h>
#include <omp.h>
#include "eigen3/Eigen/Dense"

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

//argv will accept: System (string), Nh (int), Nb (int), rs/rho (double), precision (double), degree of triples (CCD, CCDT-1, CCDT-2, etc), #threads per MPI world
int main(int argc, char** argv)
{
    Eigen::initParallel();

    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); //sets world_size equal to world size
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); //sets world_rank equal to world rank

    double       eps     = atof(argv[5]);              //remember to adjust setprecision in master when changing this
    double       conFac  = 1;                          //convergence factor
    const double pi      = M_PI;


    bool    intermediates = true;                   //turn on/off intermediates in CCD eqs
    //bool    CCDT          = true;                   //turn on/off CCDT-1
    bool    timer         = true;                   //turn on/off timer
    bool    relaxation    = true;                   //turn on/off relaxation when updating amplitudes
    double  alpha         = 0.8;                    //relaxation parameter

    bool    makeData      = false;                  //choose to write to file for a range of shells

    bool    threadsOn     = false;
    if (atoi(argv[7]) > 1){
        threadsOn = true;

        Eigen::initParallel();
        int threads = atoi(argv[7]);
        omp_set_dynamic(0);
        omp_set_num_threads(threads);
        Eigen::setNbThreads(threads);
        int n = Eigen::nbThreads( );
        std::cout << "OMP threads for world " << world_rank << ": " << n << std::endl;
    }
    else{
        omp_set_dynamic(0);
        omp_set_num_threads(8);
        if (world_rank==0){
            std::cout << "Threading is turned off, meaning:" << std::endl;
            std::cout << "- Eigen won't run parallel matrix products" << std::endl;
            std::cout << "- Permutations will be performed in serial" << std::endl;
            std::cout << "- Diagrams T5b and T5c will be run in serial" << std::endl;
        }
    }

    if (makeData == false){

        if (argc != 8 && argc != 1){
            MPI_Finalize();
            if (world_rank == 0){
                std::cout << "Failure: Missmatch command line arguments" << std::endl;
                std::cout << "The command line arguments are: Model, #particles, #shells, rs/rho, desired precision (down to 1e-16), which CC approximation" << std::endl;
                std::cout << "The CC models are: 0 -> full CCD" << std::endl;
                std::cout << "                   1 -> CCDT-1" << std::endl;
                std::cout << "                   2 -> CCDT-2" << std::endl;
                std::cout << "                   3 -> full CCDT" << std::endl;
                std::cout << "If you require any other set of CC diagrams, edit them manually in master.cpp (in function Iterator)" << std::endl;
            }
            return 0;
        }

        //we use natural units
        int Nh; int Nb;
        if (argc==8){//user defined size
            Nh = atoi(argv[2]);				//number of particles
            Nb = atoi(argv[3]);				//number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 is min for N=14)
        }
        else{        //default size
            Nh = 14;                        //number of particles
            Nb = 5;                         //number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 is min for N=14)
        }
        double  rs;     //Wigner Seitz radius
        double  rho;    //Density
        double  L3;     //Box volume
        double  L2;
        double  L1;
        double  m;      //Particle mass

        Master* master = new Master;
        master->setSize(Nh, Nb);

        //set system and physical parameters
        if (argc != 8){
            if (world_rank == 0){
                std::cout << "No arguments given, running default setup" << std::endl;
                std::cout << "Default setup: HEG for Nh=14, Nb=3, rs=1.0, 1e-16 precision, CCDT" << std::endl;
            }
            m          = 1;             //Electron mass [MeV?]
            rs         = 2;
            double  rb = 1.;            //Bohr radius [MeV^-1]
            double  r1 = pow(rs*rb, 3);
            L3         = 4.*pi*Nh*r1/3.;
            L2         = pow(L3, 2./3.);
            L1         = pow(L3, 1./3.);
            master->setSystem(new HEG(master, m, L3, L2, L1));
        }
        else if (std::string(argv[1]) == "HEG"){
            m          = 1;             //Electron mass [MeV?]
            rs         = atof(argv[4]);
            double  rb = 1.;            //Bohr radius [MeV^-1]
            double  r1 = pow(rs*rb, 3);
            L3         = 4.*pi*Nh*r1/3.;
            L2         = pow(L3, 2./3.);
            L1         = pow(L3, 1./3.);
            master->setSystem(new HEG(master, m, L3, L2, L1));
        }
        else if (std::string(argv[1]) == "MP"){
            m   = 939.565;              //Neutron mass [MeV]
            rho = atof(argv[4]);
            L3  = Nh/rho;
            L2  = pow(L3, 2./3.);
            L1  = pow(L3, 1./3.);
            master->setSystem(new MP(master, m, L3, L2, L1));
        }
        else{
            MPI_Finalize();
            if (world_rank == 0){
                std::cout << "Failure: You need to submit a model, e.g. HEG or MP" << std::endl;
                std::cout << "Which model is specified in the first command line argument" << std::endl;
            }
            return 0;
        }

        bool CCDT;

        master->setIntermediates(intermediates);
        master->setRelaxation(relaxation, alpha);
        master->setTimer(timer);
        master->setThreads(threadsOn, atoi(argv[7]));
        if (argc == 8){
            master->setCCType(atoi(argv[6]));
            if (atoi(argv[6])==0){
                master->setTriples(false);
                CCDT = false;
            }
            else{
                master->setTriples(true);
                CCDT = true;
            }
        }
        else{
            master->setCCType(3);
            master->setTriples(true);
            CCDT = true;
        }

        cout << "C++ code" << endl;

        auto t1 = Clock::now();

        if (CCDT){
            double ECCDT = master->CC_master(eps, conFac);
            if (world_rank==0){cout << "Delta ECCDT: "<< ECCDT << endl;}
        }
        else{
            double ECCD = master->CC_master(eps, conFac);
            if (world_rank==0){cout << "Delta ECCD: "<< ECCD << endl;}
        }

        auto t2 = Clock::now();

        if (world_rank == 0){
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

    }
    else if (makeData == true){

        //we use natural units

        int lower_bound = 2; int upper_bound = 10;  //set lower and upper limits on shells to be used

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
        bool CCDT = true;

        for (int j=1; j<20; j++){
            double ECC;
            double  rho      =   j/(double)100;                          //Wigner Seitz radius
            //if (world_rank==0){
                ofstream myfile;
                ostringstream os;
                os << "Nh_14_MP_rho" << rho*100 <<".txt";
                string s = os.str();
                myfile.open(s);
            //}
        for (int i=lower_bound; i<upper_bound; i++){
            int Nb = i;
            double  rb      =   1;
            //std::cout << rs << std::endl;
            //double  m       =   1;                //electron mass [MeV] (1 for HEG, 939.565 for MP)
            double  m       =   939.565;
            //double  rho     =   0.2;
            //double  r1      =   pow(rs*rb, 3);
            //double  L3      =   4.*pi*Nh*r1/3.;               //box volume
            double  L3      =   Nh/rho;
            double  L2      =   pow(L3, 2./3.);
            double  L1      =   pow(L3, 1./3.);

            Master* master = new Master;
            master->setSize(Nh, Nb);

            master->setSystem(new MP(master, m, L3, L2, L1));

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

            //if (world_rank==0){
            myfile <<  Nb << " " << "&" << " " << master->m_Ns << " " << "&" << " " << std::setprecision(16) << ECC << " \\\\ " << "\n" << std::flush;
            //}
            }

        }
        //if(world_rank==0){
        myfile.close();
        //}
    }
    }

    MPI_Finalize();

    return 0;

}