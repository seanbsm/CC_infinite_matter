#ifndef MASTER_H
#define MASTER_H

#include <Systems/system.h>
#include <Systems/heg.h>
#include <Systems/mp.h>
#include <makeampmat.h>
#include <makeintmat.h>
#include <diagrams.h>

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
    bool              m_triplesOn;
    bool              m_intermediatesOn;
    bool              m_timerOn;



    void setSize();
    void setSystem(class System* system);
    void setSize(int Nh, int Nb);
    void setIntermediates(bool argument);
    void setTriples(bool argument);
    void setTimer(bool argument);
    void setClasses();

    double CC_master(double eps, double conFac);
    double Iterator(double eps, double conFac, double E_MBPT2);
};

#endif // MASTER_H
