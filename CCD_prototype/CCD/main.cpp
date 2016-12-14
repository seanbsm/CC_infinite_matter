//imported files
//other libraries
#include <iostream>
#include <math.h>
#include <eigen3/Eigen/Dense>

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
    double  pi      =   M_PI;
    int     Nh      =   14;							//number of particles
    int     Nb      =   5;							//number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 can at max have N=14)
    double  rs      =   1;                        //Wigner Seitz radius
    double  rb      =   1;                          //Bohr radius [MeV^-1]
    double  m       =   1;                          //electron mass [MeV]
    double  L3      =   4*pi*Nh*rs/3;               //box volume
    double  L2      =   pow(L3, 2./3);
    double  L1      =   pow(L3, 1./3);

    double  eps     =   1e-10;
    double  conFac  =   1;                          //convergence factor

    Master* master = new Master;
    master->setSize(Nh, Nb);

    master->setSystem(new MP(master, m, L3, L2, L1));

    double ECCD = master->Iterator(eps, conFac);

    cout << ECCD << endl;

    return 0;
}

