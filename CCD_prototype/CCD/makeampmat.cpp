#include "makeampmat.h"
#include <iostream>

MakeAmpMat::MakeAmpMat()
{
}

void MakeAmpMat::setIntClass(class MakeIntMat* intClass){
    m_intClass = intClass;
}

void MakeAmpMat::setSystem(class System* system){
    m_system = system;
}

//recall that this function actually also needs to divide by the single fock energies (or FockMat)
void MakeAmpMat::makeBlockMat(){
    int range = m_intClass->Vhhpp.size();
    for (int i = 0; i<range; i++){
        Amplitudes.push_back( m_intClass->Vhhpp[i] );
    }
}

//for now we store the Fock matrix, but perhaps later it would be better to calculate it "on the fly"
//I deliberately declared a lot of ints here, but merely for easier debugging and readability
void MakeAmpMat::makeDenomMat(){
    for (int i=0; i<m_intClass->Vhhpp_i.size(); i++){   //remember, Vhhpp_i holds all kUnique for Vhhpp
        int lowBound_hh  = m_intClass->boundsHolder_hhpp_hh(0,i);
        int highBound_hh = m_intClass->boundsHolder_hhpp_hh(1,i);
        int lowBound_pp  = m_intClass->boundsHolder_hhpp_pp(0,i);
        int highBound_pp = m_intClass->boundsHolder_hhpp_pp(1,i);

        int dim_hh = highBound_hh - lowBound_hh;
        int dim_pp = highBound_pp - lowBound_pp;
        Eigen::MatrixXf newMat;
        newMat.conservativeResize(dim_hh, dim_pp);

        for (int hh=lowBound_hh; hh<highBound_hh; hh++){
            for (int pp=lowBound_pp; pp<highBound_pp; pp++){
                int ii = m_intClass->blockArrays_hh(1,hh);
                int jj = m_intClass->blockArrays_hh(2,hh);
                int aa = m_intClass->blockArrays_pp(1,pp);
                int bb = m_intClass->blockArrays_pp(2,pp);
                newMat(hh-lowBound_hh, pp-lowBound_pp) = 1/( (double)(m_system->f(ii) + m_system->f(jj) - m_system->f(aa) - m_system->f(bb)) );
            }
        }
        denomMat.push_back( newMat );
    }
}
