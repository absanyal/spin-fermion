#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "Observables.h"
#include "tensor_type.h"

#ifndef HF_SC_ENGINE
#define HF_SC_ENGINE


class HF_SC_Engine{
public:
    HF_SC_Engine(Parameters& Parameters__, Coordinates& Coordinates__,
             MFParams& MFParams__, Hamiltonian& Hamiltonian__,
             Observables& Observables__)
        : Parameters_(Parameters__),Coordinates_(Coordinates__),
          MFParams_(MFParams__), Hamiltonian_(Hamiltonian__),
          Observables_(Observables__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ns_(Parameters_.ns),
          ED_(Parameters_.ED_)
    {

    }

    void RUN_HF_SC_Engine();
    double Prob (double muu, double mu_new);
    double ProbCluster (double muu, double mu_new);
    Parameters &Parameters_;
    Coordinates &Coordinates_;
    MFParams &MFParams_;
    Hamiltonian &Hamiltonian_;
    Observables &Observables_;
    const int lx_, ly_, ns_;
    bool ED_;

};

/*
 * ***********
 *  Functions in Class MCEngine ------
 *  ***********
*/

void HF_SC_Engine::RUN_HF_SC_Engine(){

    complex<double> zero(0.0,0.0);
    bool Simple_Mixing = Parameters_.Simple_Mixing;








} // ---------



#endif // HF_SC_ENGINE
