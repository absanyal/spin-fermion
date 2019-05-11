#ifndef Parameters_class
#define Parameters_class
#include "tensor_type.h"

class Parameters{

public:
    int lx, ly, ns, IterMax, MCNorm, RandomSeed;
    int TBC_mx, TBC_my;
    int TBC_cellsX, TBC_cellsY;
    int lx_cluster, ly_cluster;
    double mus, mus_Cluster, Fill,pi;
    double J_Hund, K1x ,K1y;
    double k_const;
    double Disorder_Strength, RandomDisorderSeed;

    double t_hopping;

    Mat_1_string MC_DOF;


    bool Cooling_;
    bool ED_;

    bool Metropolis_Algorithm;
    bool Heat_Bath_Algorithm;
    bool Perform_HF_SC_calculation;
    bool Simple_Mixing;

    bool MC_on_theta, MC_on_phi, MC_on_moment_size, MC_on_local_density;

    /*
SavingMicroscopicStates=1
NoOfMicroscopicStates=50
      */
    bool Saving_Microscopic_States;
    int No_Of_Microscopic_States;

    double temp_max, beta_min;
    double temp_min, beta_max;
    double d_Temp;

    int Last_n_sweeps_for_measurement;
    int Measurement_after_each_m_sweeps;

    double temp,beta,Eav,maxmoment;
    double WindowSize, AccCount[2];
    char Dflag;

    void Initialize(string inputfile_);
    double matchstring(string file,string match);
    string matchstring2(string file,string match);

};


void Parameters::Initialize(string inputfile_){

    maxmoment=1.0;
    t_hopping=1.0;
    string monte_carlo_variables_;
    double cooling_double;
    double metropolis_double;
    double Perform_HF_SC_calculation_double;
    double Simple_Mixing_double;
    double ED_double;
    int SavingMicroscopicStates_int;
    string temp_string;


    cout << "____________________________________" << endl;
    cout << "Reading the inputfile: " << inputfile_ << endl;
    cout << "____________________________________" << endl;

    monte_carlo_variables_=matchstring2(inputfile_,"Monte_Carlo_variables");

    stringstream monte_carlo_variables_stream(monte_carlo_variables_);
    int no_of_MC_vars;
    monte_carlo_variables_stream >> no_of_MC_vars;
    assert(no_of_MC_vars<5);
    MC_on_theta=false;
    MC_on_phi=false;
    MC_on_moment_size=false;
    MC_on_local_density=false;
    MC_DOF.resize(no_of_MC_vars);

    for(int i=0;i<no_of_MC_vars;i++){
        monte_carlo_variables_stream >> temp_string;
        MC_DOF[i]=temp_string;

        if(temp_string=="phi"){
            MC_on_phi=true;
        }
        else if(temp_string=="theta"){
            MC_on_theta=true;
        }
        else if(temp_string=="moment_size"){
            MC_on_moment_size=true;
        }
        else if(temp_string=="local_density"){
            MC_on_local_density=true;
        }
        else{
            cout<<"MC cannot be performed on "<<temp_string<<
                  ", please choose only from phi theta moment_size local_density"<<endl;
            assert(false);
        }
    }




    lx = int(matchstring(inputfile_,"Xsite"));
    ly = int(matchstring(inputfile_,"Ysite"));
    TBC_mx = int(matchstring(inputfile_,"TwistedBoundaryCond_mx"));
    TBC_my = int(matchstring(inputfile_,"TwistedBoundaryCond_my"));
    TBC_cellsX = int(matchstring(inputfile_,"TBC_cellsX"));
    TBC_cellsY = int(matchstring(inputfile_,"TBC_cellsY"));
    lx_cluster = int(matchstring(inputfile_,"Cluster_lx"));
    ly_cluster = int(matchstring(inputfile_,"Cluster_ly"));
    SavingMicroscopicStates_int = int (matchstring(inputfile_,"SavingMicroscopicStates"));

    assert(SavingMicroscopicStates_int==1 ||
           SavingMicroscopicStates_int==0);
    if(SavingMicroscopicStates_int==1){
        Saving_Microscopic_States=true;
    }
    else{
        Saving_Microscopic_States=false;
    }


    No_Of_Microscopic_States = int(matchstring(inputfile_,"NoOfMicroscopicStates"));



    ns = lx*ly;
    cout << "TotalNumberOfSites = "<< ns << endl;

    Fill = matchstring(inputfile_,"Fill");
    cout << "TotalNumberOfParticles = "<< ns*Fill*2.0 << endl;

    IterMax = int(matchstring(inputfile_,"MaxMCsweeps"));
    MCNorm = 0.0; //matchstring(inputfile,"MCNorm")
    RandomSeed = matchstring(inputfile_,"RandomSeed");
    RandomDisorderSeed = matchstring(inputfile_,"RandomDisorderSeed");
    Disorder_Strength = matchstring(inputfile_,"Disorder_Strength");
    J_Hund = matchstring(inputfile_,"J_Hund");
    K1x = matchstring(inputfile_,"K");
    K1y=K1x;

    Dflag = 'N';

    metropolis_double=double(matchstring(inputfile_,"Metropolis_Algo"));
    if(metropolis_double==1.0){
        Metropolis_Algorithm=true;
        Heat_Bath_Algorithm=false;
    }
    else if(metropolis_double==0.0){
        Metropolis_Algorithm=false;
        Heat_Bath_Algorithm=true;
    }
    else{
        cout<<"ERROR: Metropolis_Algo can be only 1 (true) or 0 (false)"<<endl;
        assert(metropolis_double==0.0);
    }

    Perform_HF_SC_calculation_double=double(matchstring(inputfile_,"Perform_HF_SC_calculation"));
    if(Perform_HF_SC_calculation_double==1.0){
        Perform_HF_SC_calculation=true;
    }
    else{
        Perform_HF_SC_calculation=false;
    }

    Simple_Mixing_double==double(matchstring(inputfile_,"Simple_Mixing"));
    if(Simple_Mixing_double==1.0){
        Simple_Mixing=true;
    }
    else{
        Simple_Mixing=false;
    }


    cooling_double=double(matchstring(inputfile_,"Cooling"));
    if(cooling_double==1.0){
        Cooling_=true;

        temp_min = double(matchstring(inputfile_,"Temperature_min"));
        temp_max = double(matchstring(inputfile_,"Temperature_max"));
        d_Temp = double(matchstring(inputfile_,"dTemperature"));
        beta_max=double(11604.0/ temp_min);
        beta_min=double(11604.0/ temp_max);

    }
    else if(cooling_double==0.0){
        Cooling_=false;

        temp = double(matchstring(inputfile_,"Temperature"));   // temperature in kelvin
        beta=double(11604.0/ temp);    //Beta which is (T*k_b)^-1

        temp_min = temp;
        temp_max = temp;
        d_Temp=10.0;//arbitrary positive number

    }
    else{
        cout<<"ERROR: Cooling can be only 1 (true) or 0 (false)"<<endl;
        assert(cooling_double==0.0);
    }



    ED_double=double(matchstring(inputfile_,"Perform_ED"));
    if(ED_double==1.0){
        ED_=true;
        lx_cluster=lx;
        ly_cluster=ly;
    }
    else if(ED_double==0.0){
        ED_=false;
    }
    else{
        cout<<"ERROR: Perform_ED can be only 1 (true) or 0 (false)"<<endl;
        assert(ED_double==0.0);
    }



    Last_n_sweeps_for_measurement=int(matchstring(inputfile_,"Last_n_sweeps_for_measurement"));
    Measurement_after_each_m_sweeps=int(matchstring(inputfile_,"Measurement_after_each_m_sweeps"));


    pi=4.00*atan(double(1.0));
    Eav=0.0;
    AccCount[0]=0;
    AccCount[1]=0;

    WindowSize=double(0.01);
    mus=0.25;
    cout << "____________________________________" << endl;
}


double Parameters::matchstring(string file,string match) {
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass=false;
    while (std::getline(readFile, line)) {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass==false) {
            // ---------------------------------
            if (iss >> amount && test==match) {
                // cout << amount << endl;
                pass=true;
            }
            else {
                pass=false;
            }
            // ---------------------------------
            if(pass) break;
        }
    }
    if (pass==false) {
        string errorout=match;
        errorout+="= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

string Parameters::matchstring2(string file,string match) {

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    if(readFile.is_open())
    {
        while(!readFile.eof())
        {
            getline(readFile,line);

            if ((offset = line.find(match, 0)) != string::npos) {
                amount = line.substr (offset+match.length()+1);				}

        }
        readFile.close();
    }
    else
    {cout<<"Unable to open input file while in the Parameters class."<<endl;}




    cout << match << " = " << amount << endl;
    return amount;
}

#endif



