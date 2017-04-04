#ifndef MASTER_H
#define MASTER_H

#include <Systems/system.h>


class Master
{
public:
    class System*     m_system             = nullptr;
    class MakeAmpMat* m_ampClass           = nullptr;
    class MakeIntMat* m_intClass           = nullptr;
    class Diagrams*   m_diagrams           = nullptr;
    int               m_Nh                 = 0;        //number of particles
    int               m_Nb                 = 0;
    int               m_Ns                 = 0;
    int               m_Np                 = 0;
    bool              m_triplesOn          = false;
    bool              m_intermediatesOn    = false;
    bool              m_timerOn            = false;
    bool              m_relaxation         = false;
    double            m_alpha              = 1;

    int world_rank;

    void setSize();
    void setSystem(class System* system);
    void setSize(int Nh, int Nb);
    void setIntermediates(bool argument);
    void setTriples(bool argument);
    void setTimer(bool argument);
    void setRelaxation(bool argument, double alpha);
    void setClasses();

    double CC_master(double eps, double conFac);
    double Iterator(double eps, double conFac, double E_MBPT2);
};

#endif // MASTER_H
