#include "makeampmat.h"
#include <iostream>

MakeAmpMat::MakeAmpMat()
{
}

//recall that this function actually also needs to divide by the single fock energies (or FockMat)
void MakeAmpMat::makeBlockMat(std::vector<Eigen::MatrixXf> Vhhpp){
    int range = Vhhpp.size();
    for (int i = 0; i<range; i++){
        Amplitudes.push_back( Vhhpp[i] );
    }
}
