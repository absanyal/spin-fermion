#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "random"
#include <stdlib.h>
#define PI acos(-1.0)

#ifndef MFParams_class
#define MFParams_class

class MFParams{ 
public:

    // Define Fields
    Matrix<double> etheta, ephi;
    Matrix<double> Sz, Sx, Sy;
    Matrix<double> etheta_avg, ephi_avg;
    Matrix<double> Moment_Size;
    Matrix<double> Local_density;
    Matrix<double> Disorder;

    // Constructor
    MFParams(Parameters& Parameters__, Coordinates&  Coordinates__, mt19937_64& Generator1__ , mt19937_64& Generator2__)
        :Parameters_(Parameters__),Coordinates_(Coordinates__), Generator1_(Generator1__), Generator2_(Generator2__)
    {
        //setupm_arr();
        initialize();
    }


    double random1();
    double random2();
    void FieldThrow(int site, string mc_dof_type);
    void initialize();
    void Adjust_MCWindow();
    void Calculate_Fields_Avg();
    void Read_classical_DOFs(string filename);


    Parameters &Parameters_;
    Coordinates &Coordinates_;
    mt19937_64 &Generator1_; //for random fields
    mt19937_64 &Generator2_; //for random disorder
    int lx_,ly_,ns_;

    uniform_real_distribution<double> dis1_;//for random fields
    uniform_real_distribution<double> dis2_;//for random disorder

    //mt19937_64 mt_rand(Parameters_.RandomSeed);


};


void MFParams::Adjust_MCWindow(){
    double ratio;
    ratio=Parameters_.AccCount[0]/(Parameters_.AccCount[0]+Parameters_.AccCount[1]);
    //cout<<"ratio= "<< ratio << "temp= "<<Parameters_.temp << endl;
    Parameters_.AccCount[0]=0; Parameters_.AccCount[1]=0;
    Parameters_.WindowSize *= abs(1.0 + 1.0*(ratio-0.5));
    //Parameters_.WindowSize =1.0;
    cout << "Ratio: " << ratio << "  window size:  "<<Parameters_.WindowSize<< endl;
    return;
} // ----------


void MFParams::FieldThrow(int site, string mc_dof_type){
    int a,b;

    double Pi=Parameters_.pi;
    double MC_Window = Parameters_.WindowSize;

    a=Coordinates_.indx(site);
    b=Coordinates_.indy(site);


    //ANGLES
    if(mc_dof_type=="phi"){
        ephi(a,b) += 2*Pi*(random1()-0.5)*MC_Window;
        if( ephi(a,b) < 0.0) {ephi(a,b) += 2.0*Pi; }
        if( ephi(a,b) >=2.0*Pi) {ephi(a,b) -= 2.0*Pi;}
    }

    if(mc_dof_type=="theta"){
        etheta(a,b) +=  Pi*(random1()-0.5)*MC_Window;
        if ( etheta(a,b) < 0.0 ) {
            etheta(a,b) = - etheta(a,b);
            ephi(a,b) = fmod( ephi(a,b)+Pi, 2.0*Pi );
        }

        if ( etheta(a,b) > Pi ) {
            etheta(a,b) -= 2.0*Pi ;
            ephi(a,b) = fmod( ephi(a,b) + Pi, 2.0*Pi );
        }

    }

    //Moment Size
    if(mc_dof_type=="moment_size"){
        if(Parameters_.MC_on_moment_size == true){
            Moment_Size(a,b) = abs(Moment_Size(a,b) + ((random1()-0.5)*MC_Window));
        }
    }


    // Local Density
    if(mc_dof_type=="local_density"){
        if(Parameters_.MC_on_local_density == true){
            Local_density(a,b) = abs( Local_density(a,b)  +  ((random1()-0.5)*MC_Window)    );
        }
    }



} // ----------


double MFParams::random1(){

    return dis1_(Generator1_);

}

double MFParams::random2(){

    return dis2_(Generator2_);

}


void MFParams::initialize(){

    lx_=Coordinates_.lx_;
    ly_=Coordinates_.ly_;

    // srand(Parameters_.RandomSeed);


    etheta_avg.resize(lx_,ly_);
    ephi_avg.resize(lx_,ly_);
    Disorder.resize(lx_,ly_);

    etheta.resize(lx_,ly_);
    ephi.resize(lx_,ly_);
    Moment_Size.resize(lx_,ly_);
    Local_density.resize(lx_,ly_);

    Sz.resize(lx_,ly_);
    Sx.resize(lx_,ly_);
    Sy.resize(lx_,ly_);

    ofstream Disorder_conf_file("Disorder_conf_used");
    Disorder_conf_file<<"#seed="<<Parameters_.RandomDisorderSeed<<
                        " for mt19937_64 Generator is used"<<endl;
    Disorder_conf_file<<"#ix   iy    Dis[ix,iy]"<<endl;

    ofstream Initial_MC_DOF_file("Initial_MC_DOF_values");
    ofstream Initial_HF_OrderParameters_file("Initial_HF_OrderParameters_values");


    if(Parameters_.Perform_HF_SC_calculation==false){

        Initial_MC_DOF_file<<"#seed="<<Parameters_.RandomSeed<<
                             " for mt19937_64 Generator is used"<<endl;
        Initial_MC_DOF_file<<"#ix   iy   Theta(x,y)    Phi(x,y)       Moment_Size(x,y)      Local_density(x,y)"<<endl;
    }
    else if(Parameters_.Perform_HF_SC_calculation==true){

        Initial_HF_OrderParameters_file<<"#seed="<<Parameters_.RandomSeed<<
                                         " for mt19937_64 Generator is used"<<endl;
        Initial_HF_OrderParameters_file<<"#ix   iy   Sz(x,y)    Sx(x,y)   Sy(x,y)   Local_density(x,y)"<<endl;
    }



    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){
            //ephi(i,j)=(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
            //etheta(i,j)=0.5*Parameters_.pi + grnd()*0.2;

            //q=(pi,pi)
            // ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
            // etheta(i,j)=0.5*(pow(-1.0,j+i)  + 1.0 )*PI ;//+ grnd()*0.2;

            //q=(0,pi)
            //ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
            //etheta(i,j)=0.5*(pow(-1.0,j)  + 1.0 )*PI; //+ grnd()*0.2;

            //q=(0,0)
            // ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
            // etheta(i,j)=0.0; //+ grnd()*0.2;


            if(Parameters_.Perform_HF_SC_calculation==false){
                //RANDOM fields
                if(Parameters_.MC_on_phi == true){
                    ephi(i,j)=2.0*random1()*PI;
                }
                else{
                    ephi(i,j)=0.0;
                }

                if(Parameters_.MC_on_theta == true){
                    etheta(i,j)=random1()*PI;
                }
                else{
                    etheta(i,j)=0.0;
                }

                if(Parameters_.MC_on_moment_size == true){
                    Moment_Size(i,j) =  1.0;//random1();
                    //   Moment_Size(i,j)=0.5;
                    //Moment_Size(i,j)=1.0 - (random1() - 0.5)*0.1;
                }
                else{
                    Moment_Size(i,j)=1.0;
                }

                if(Parameters_.MC_on_local_density == true){
                    Local_density(i,j)=random1();
                }
                else{
                    Local_density(i,j)=Parameters_.Fill*2.0;
                }
            }
            else if(Parameters_.Perform_HF_SC_calculation==true){
                Sz(i,j)=random1();
                Sx(i,j)=random1();
                Sy(i,j)=random1();
                Local_density(i,j)=random1();
            }




            //RANDOM Disorder
            Disorder(i,j)=Parameters_.Disorder_Strength*((2.0*random2())-1.0);
            Disorder_conf_file<<i<<"  "<<j<<"  "<<Disorder(i,j)<<endl;


            if(Parameters_.Perform_HF_SC_calculation==false){
            Initial_MC_DOF_file<<i<<setw(15)<<j<<setw(15)<<etheta(i,j)<<setw(15)<<ephi(i,j)
                              <<setw(15)<<Moment_Size(i,j)<<setw(15)<<Local_density(i,j)<<endl;
            }
            else if(Parameters_.Perform_HF_SC_calculation==true){
            Initial_HF_OrderParameters_file<<i<<setw(15)<<j<<setw(15)<<Sz(i,j)<<setw(15)<<Sx(i,j)
                                          <<setw(15)<<Sy(i,j)<<setw(15)<<Local_density(i,j)<<endl;
            }



        }
        Disorder_conf_file<<endl;
        Initial_HF_OrderParameters_file<<endl;

    }


} // ----------

void MFParams::Calculate_Fields_Avg(){

    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){

            ephi_avg(i,j)= ephi_avg(i,j) + ephi(i,j);
            etheta_avg(i,j)=etheta_avg(i,j) + etheta(i,j) ;
        }
    }


} // ----------


void MFParams::Read_classical_DOFs(string filename){

    string tmp_str;
    double tmp_double;
    ifstream fl_in(filename.c_str());
    fl_in >> tmp_str;

    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            fl_in>>tmp_double>>tmp_double>>etheta(i,j)>>ephi(i,j);
        }
    }


} // ----------


#endif
