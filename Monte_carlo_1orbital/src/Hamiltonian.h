#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#define PI acos(-1.0)

#ifndef Hamiltonian_class
#define Hamiltonian_class

extern "C" void zheev_(char *, char *, int *, std::complex<double> *, int *, double *,
                       std::complex<double> *, int *, double *, int *);

class Hamiltonian
{
public:
    Hamiltonian(Parameters &Parameters__, Coordinates &Coordinates__, Coordinates &CoordinatesCluster__, MFParams &MFParams__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), CoordinatesCluster_(CoordinatesCluster__), MFParams_(MFParams__)

    {
        Initialize();
        Hoppings();
        HTBCreate();
        HTBClusterCreate();
    }

    void Initialize();                                     //::DONE
    void Hoppings();                                       //::DONE
    double GetCLEnergy();                                  //::DONE
    void InteractionsCreate();                             //::DONE
    void InteractionsClusterCreate(int Center_site);       //::DONE
    void Check_Hermiticity();                              //::DONE
    void Check_up_down_symmetry();                         //::DONE
    void HTBCreate();                                      //::DONE
    void HTBClusterCreate();                               //::DONE
    double chemicalpotential(double muin, double filling); //::DONE

    double chemicalpotentialCluster(double muin, double filling); //::DONE

    double TotalDensity();   //::DONE
    double ClusterDensity(); //::DONE
    double E_QM();           //::DONE

    double E_QMCluster();                 //::DONE
    void Diagonalize(char option);        //::DONE
    void DiagonalizeCluster(char option); //::DONE
    void copy_eigs(int i);                //::DONE
    void copy_eigs_Cluster(int i);        //::DONE

    Parameters &Parameters_;
    Coordinates &Coordinates_;
    Coordinates &CoordinatesCluster_;
    MFParams &MFParams_;
    int lx_, ly_, ns_, orbs_;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> HTBCluster_;
    Matrix<complex<double>> Ham_;
    Matrix<complex<double>> HamCluster_;
    Matrix<double> Tx, Ty, Tpxpy, Tpxmy;
    vector<double> eigs_, eigsCluster_, eigsCluster_saved_, eigs_saved_, sx_, sy_, sz_;

    double HS_factor;
};

double Hamiltonian::chemicalpotential(double muin, double filling)
{
    double mu_out;
    double n1, N;
    double dMubydN;
    double nstate = eigs_.size();
    dMubydN = 0.05 * (eigs_[nstate - 1] - eigs_[0]) / nstate;
    N = filling * double(eigs_.size());
    //temp=Parameters_.temp;
    mu_out = muin;
    bool converged = false;

    if (1 == 2)
    {
        for (int i = 0; i < 100000; i++)
        {
            n1 = 0.0;
            for (int j = 0; j < nstate; j++)
            {
                n1 += double(1.0 / (exp((eigs_[j] - mu_out) * Parameters_.beta) + 1.0));
            }
            //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
            if (abs(N - n1) < double(0.0001))
            {
                //cout<<abs(N-n1)<<endl;
                converged = true;
                break;
            }
            else
            {
                mu_out += (N - n1) * dMubydN;
                //cout<<i<<"    "<<n1<<"    "<<N-n1<<endl;
            }
        }

        if (!converged)
        {
            cout<<"mu_not_converged, N = "<<n1<<endl;
        }
        else
        {
            //cout<<"mu converged, N = "<<n1<<endl;
        }
    }

    double mu1, mu2;
    double mu_temp = muin;
    //cout<<"mu_input = "<<mu_temp<<endl;
    if (1 == 1)
    {
        mu1 = eigs_[0];
        mu2 = eigs_[nstate - 1];
        for (int i = 0; i < 40000; i++)
        {
            n1 = 0.0;
            for (int j = 0; j < nstate; j++)
            {
                n1 += double(1.0 / (exp((eigs_[j] - mu_temp) * Parameters_.beta) + 1.0));
            }
            //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
            if (abs(N - n1) < double(0.0001))
            {
                //cout<<abs(N-n1)<<endl;
                converged = true;
                break;
            }
            else
            {
                if (n1 > N)
                {
                    mu2 = mu_temp;
                    mu_temp = 0.5 * (mu1 + mu_temp);
                }
                else
                {
                    mu1 = mu_temp;
                    mu_temp = 0.5 * (mu2 + mu_temp);
                }
            }
            //cout<<"mu_temp = "<<mu_temp<<"   "<<n1<<endl;
        }

        if (!converged)
        {
            cout<<"mu_not_converged, N = "<<n1<<endl;
        }
        else
        {
            //cout<<"mu converged, N = "<<n1<<endl;
        }

        mu_out = mu_temp;
    }

    return mu_out;
} // ----------

double Hamiltonian::chemicalpotentialCluster(double muin, double filling)
{
    double mu_out;
    double n1, N;
    double dMubydN;
    double nstate = eigsCluster_.size();
    dMubydN = 0.05 * (eigsCluster_[nstate - 1] - eigsCluster_[0]) / nstate;
    N = filling * double(eigsCluster_.size());
    //temp=Parameters_.temp;
    mu_out = muin;
    bool converged = false;

    if (1 == 2)
    {
        for (int i = 0; i < 100000; i++)
        {
            n1 = 0.0;
            for (int j = 0; j < nstate; j++)
            {
                n1 += double(1.0 / (exp((eigsCluster_[j] - mu_out) * Parameters_.beta) + 1.0));
            }
            //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
            if (abs(N - n1) < double(0.0001))
            {
                //cout<<abs(N-n1)<<endl;
                converged = true;
                break;
            }
            else
            {
                mu_out += (N - n1) * dMubydN;
                //cout<<i<<"    "<<n1<<"    "<<N-n1<<endl;
            }
        }

        if (!converged)
        {
            cout<<"mu_not_converged, N = "<<n1<<endl;
        }
        else
        {
            //cout<<"mu converged, N = "<<n1<<endl;
        }
    }

    double mu1, mu2;
    double mu_temp = muin;
    //cout<<"mu_input = "<<mu_temp<<endl;
    if (1 == 1)
    {
        mu1 = eigsCluster_[0];
        mu2 = eigsCluster_[nstate - 1];
        for (int i = 0; i < 40000; i++)
        {
            n1 = 0.0;
            for (int j = 0; j < nstate; j++)
            {
                n1 += double(1.0 / (exp((eigsCluster_[j] - mu_temp) * Parameters_.beta) + 1.0));
            }
            //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
            if (abs(N - n1) < double(0.0001))
            {
                //cout<<abs(N-n1)<<endl;
                converged = true;
                break;
            }
            else
            {
                if (n1 > N)
                {
                    mu2 = mu_temp;
                    mu_temp = 0.5 * (mu1 + mu_temp);
                }
                else
                {
                    mu1 = mu_temp;
                    mu_temp = 0.5 * (mu2 + mu_temp);
                }
            }
            //cout<<"mu_temp = "<<mu_temp<<"   "<<n1<<endl;
        }

        if (!converged)
        {
            cout<<"mu_not_converged, N = "<<n1<<endl;
        }
        else
        {
            //cout<<"mu converged, N = "<<n1<<endl;
        }

        mu_out = mu_temp;
    }

    return mu_out;
} // ----------

void Hamiltonian::Initialize()
{

    //For Hubbard Stratonovich transformation
    HS_factor = 0.0;

    //else use
    //HS_factor=0.0;

    int ns = (Parameters_.lx_cluster) * (Parameters_.ly_cluster);

    ly_ = Parameters_.ly;
    lx_ = Parameters_.lx;
    ns_ = Parameters_.ns;

    int space = 2 * ns_;
    int spaceCluster = 2 * ns;

    HTB_.resize(space, space);
    Ham_.resize(space, space);
    HTBCluster_.resize(spaceCluster, spaceCluster);
    HamCluster_.resize(spaceCluster, spaceCluster);
    eigs_.resize(space);
    sx_.resize(space);
    sy_.resize(space);
    sz_.resize(space);
    eigs_saved_.resize(space);
    eigsCluster_.resize(spaceCluster);
    eigsCluster_saved_.resize(spaceCluster);

} // ----------

double Hamiltonian::TotalDensity()
{
    double n1 = 0.0;
    for (int j = 0; j < eigs_.size(); j++)
    {
        n1 += 1.0f / (exp(Parameters_.beta * (eigs_[j] - Parameters_.mus)) + 1.0);
    }
    return n1;
} // ----------

double Hamiltonian::ClusterDensity()
{
    double n1 = 0.0;
    for (int j = 0; j < eigsCluster_.size(); j++)
    {
        n1 += 1.0f / (exp(Parameters_.beta * (eigsCluster_[j] - Parameters_.mus_Cluster)) + 1.0);
    }
    return n1;
} // ----------

double Hamiltonian::E_QM()
{
    double E = 0.0;
    for (int j = 0; j < eigs_.size(); j++)
    {
        //E +=  (eigs_[j]-Parameters_.mus)/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
        E += (eigs_[j]) / (exp(Parameters_.beta * (eigs_[j] - Parameters_.mus)) + 1.0);
    }
    return E;
} // ----------

double Hamiltonian::E_QMCluster()
{
    double E = 0.0;
    for (int j = 0; j < eigsCluster_.size(); j++)
    {
        //E +=  (eigs_[j]-Parameters_.mus)/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
        E += (eigsCluster_[j]) / (exp(Parameters_.beta * (eigsCluster_[j] - Parameters_.mus_Cluster)) + 1.0);
    }
    return E;
} // ----------

double Hamiltonian::GetCLEnergy()
{

    double EClassical;
    int site;
    double ei, ai;

    for (int i = 0; i < lx_; i++)
    {
        for (int j = 0; j < ly_; j++)
        {
            site = Coordinates_.Nc(i, j); //+x
            ei = MFParams_.etheta(i, j);
            ai = MFParams_.ephi(i, j);
            sx_[site] = MFParams_.Moment_Size(i, j) * cos(ai) * sin(ei);
            sy_[site] = MFParams_.Moment_Size(i, j) * sin(ai) * sin(ei);
            sz_[site] = MFParams_.Moment_Size(i, j) * cos(ei);
        }
    }

    // Classical Energy
    EClassical = double(0.0);

    for (int ix = 0; ix < lx_; ix++)
    {
        for (int iy = 0; iy < ly_; iy++)
        {

            EClassical += HS_factor * (-0.5) * Parameters_.J_Hund * ((MFParams_.Moment_Size(ix, iy) * MFParams_.Moment_Size(ix, iy)) - (0.25 * MFParams_.Local_density(ix, iy) * MFParams_.Local_density(ix, iy)));
        }
    }

    for (int i = 0; i < ns_; i++)
    {
        site = Coordinates_.neigh(i, 0); //+x
        EClassical += 1.0 * Parameters_.K1x * (sx_[i] * sx_[site] + sy_[i] * sy_[site] + 1.0 * sz_[i] * sz_[site]);
        site = Coordinates_.neigh(i, 2); //+y
        EClassical += Parameters_.K1y * (sx_[i] * sx_[site] + sy_[i] * sy_[site] + 1.0 * sz_[i] * sz_[site]);
    }

    return EClassical;
} // ----------

void Hamiltonian::InteractionsCreate()
{

    /*
    For Effective Hamiltonian derived from
    Hubbard-Stratonovich transformation, remember:
    JH= -2*U
    HS_factor=1.0;
    K=0;
     */

    int space = 2 * ns_;
    int a;
    double ei, ai;
    double den;

    Ham_ = HTB_;
    // Ham_.print();

    for (int i = 0; i < ns_; i++)
    { // For each site
        ei = MFParams_.etheta(Coordinates_.indx(i), Coordinates_.indy(i));
        ai = MFParams_.ephi(Coordinates_.indx(i), Coordinates_.indy(i));
        den = MFParams_.Local_density(Coordinates_.indx(i), Coordinates_.indy(i));

        Ham_(i, i) += HS_factor * (-0.25) * Parameters_.J_Hund * (den);
        Ham_(i + ns_, i + ns_) += HS_factor * (-0.25) * Parameters_.J_Hund * (den);
        Ham_(i, i) += Parameters_.J_Hund * (cos(ei)) * 0.5 * MFParams_.Moment_Size(Coordinates_.indx(i), Coordinates_.indy(i));
        Ham_(i + ns_, i + ns_) += Parameters_.J_Hund * (-cos(ei)) * 0.5 * MFParams_.Moment_Size(Coordinates_.indx(i), Coordinates_.indy(i));
        Ham_(i, i + ns_) += Parameters_.J_Hund * sin(ei) * complex<double>(cos(ai), -sin(ai)) * 0.5 * MFParams_.Moment_Size(Coordinates_.indx(i), Coordinates_.indy(i)); //S-
        Ham_(i + ns_, i) += Parameters_.J_Hund * sin(ei) * complex<double>(cos(ai), sin(ai)) * 0.5 * MFParams_.Moment_Size(Coordinates_.indx(i), Coordinates_.indy(i));  //S+

        // On-Site potential
        for (int spin = 0; spin < 2; spin++)
        {
            a = i + ns_ * spin;
            Ham_(a, a) = complex<double>(1.0, 0.0) * MFParams_.Disorder(Coordinates_.indx(i), Coordinates_.indy(i));
        }
    }

} // ----------

void Hamiltonian::InteractionsClusterCreate(int Center_site)
{

    int ns = (Parameters_.lx_cluster) * (Parameters_.ly_cluster);

    int space = 2 * ns;
    int x_pos, y_pos;
    double ei, ai, den;
    int a;

    HamCluster_ = HTBCluster_;

    for (int i = 0; i < ns; i++)
    { // For each site in cluster
        x_pos = Coordinates_.indx(Center_site) - int(Parameters_.lx_cluster / 2) + CoordinatesCluster_.indx(i);
        y_pos = Coordinates_.indy(Center_site) - int(Parameters_.ly_cluster / 2) + CoordinatesCluster_.indy(i);
        x_pos = (x_pos + Coordinates_.lx_) % Coordinates_.lx_;
        y_pos = (y_pos + Coordinates_.ly_) % Coordinates_.ly_;

        ei = MFParams_.etheta(x_pos, y_pos);
        ai = MFParams_.ephi(x_pos, y_pos);
        den = MFParams_.Local_density(x_pos, y_pos);

        HamCluster_(i, i) += HS_factor * (-0.25) * Parameters_.J_Hund * (den);
        HamCluster_(i + ns, i + ns) += HS_factor * (-0.25) * Parameters_.J_Hund * (den);
        HamCluster_(i, i) += Parameters_.J_Hund * (cos(ei)) * 0.5 * MFParams_.Moment_Size(x_pos, y_pos);
        HamCluster_(i + ns, i + ns) += Parameters_.J_Hund * (-cos(ei)) * 0.5 * MFParams_.Moment_Size(x_pos, y_pos);
        HamCluster_(i, i + ns) += Parameters_.J_Hund * sin(ei) * complex<double>(cos(ai), -sin(ai)) * 0.5 * MFParams_.Moment_Size(x_pos, y_pos); //S-
        HamCluster_(i + ns, i) += Parameters_.J_Hund * sin(ei) * complex<double>(cos(ai), sin(ai)) * 0.5 * MFParams_.Moment_Size(x_pos, y_pos);  //S+

        for (int spin = 0; spin < 2; spin++)
        {
            a = i + ns * spin;
            HamCluster_(a, a) = complex<double>(1.0, 0.0) * MFParams_.Disorder(x_pos, y_pos);
        }
    }

} // ----------

void Hamiltonian::Check_up_down_symmetry()

{
    complex<double> temp(0, 0);
    complex<double> temp2;

    for (int i = 0; i < ns_; i++)
    {
        for (int j = 0; j < ns_; j++)
        {
            temp2 = Ham_(i, j) - Ham_(i + ns_, j + ns_); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            temp += temp2 * conj(temp2);
        }
    }

    cout << "Assymetry in up-down sector: " << temp << endl;
}

void Hamiltonian::Check_Hermiticity()

{
    complex<double> temp(0, 0);
    complex<double> temp2;

    for (int i = 0; i < HamCluster_.n_row(); i++)
    {
        for (int j = 0; j < HamCluster_.n_row(); j++)
        {
            if (HamCluster_(i, j) != conj(HamCluster_(j, i)))
            {
                cout << i << "," << j << endl;
                cout << "i,j = " << HamCluster_(i, j) << ", j,i=" << conj(HamCluster_(j, i)) << endl;
            }
            assert(HamCluster_(i, j) == conj(HamCluster_(j, i))); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            //temp +=temp2*conj(temp2);
        }
    }

    // cout<<"Hermiticity: "<<temp<<endl;
}

void Hamiltonian::Diagonalize(char option)
{

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);

    char jobz = option;
    // jobz = 'V';
    char uplo = 'U'; //WHY ONLY 'L' WORKS?
    int n = Ham_.n_row();
    int lda = Ham_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3 * n - 2);
    int info;
    int lwork = -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(), eigs_.end(), 0);
    // query:
    zheev_(&jobz, &uplo, &n, &(Ham_(0, 0)), &lda, &(eigs_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz, &uplo, &n, &(Ham_(0, 0)), &lda, &(eigs_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    if (info != 0)
    {
        std::cerr << "info=" << info << "\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}
}

void Hamiltonian::DiagonalizeCluster(char option)
{

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);

    char jobz = option;
    // jobz = 'V';
    // cout << option;
    char uplo = 'U'; //WHY ONLY 'L' WORKS?
    int n = HamCluster_.n_row();
    int lda = HamCluster_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3 * n - 2);
    int info;
    int lwork = -1;

    eigsCluster_.resize(HamCluster_.n_row());
    fill(eigsCluster_.begin(), eigsCluster_.end(), 0);
    // query:
    zheev_(&jobz, &uplo, &n, &(HamCluster_(0, 0)), &lda, &(eigsCluster_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz, &uplo, &n, &(HamCluster_(0, 0)), &lda, &(eigsCluster_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    if (info != 0)
    {
        std::cerr << "info=" << info << "\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}
}

void Hamiltonian::HTBCreate()
{

    int mx = Parameters_.TBC_mx;
    int my = Parameters_.TBC_my;
    complex<double> phasex, phasey;
    int l, m, a, b;

    HTB_.fill(0.0);

    for (l = 0; l < ns_; l++)
    {

        // * +x direction Neighbor
        if (Coordinates_.indx(l) == (Coordinates_.lx_ - 1))
        {
            phasex = exp(iota_complex * 2.0 * (1.0 * mx) * PI / (1.0 * Parameters_.TBC_cellsX));
            phasey = one_complex;
        }
        else
        {
            phasex = one_complex;
            phasey = one_complex;
        }
        m = Coordinates_.neigh(l, 0);
        for (int spin = 0; spin < 2; spin++)
        {

            a = l + ns_ * spin;
            b = m + ns_ * spin;
            assert(a != b);
            if (a != b)
            {
                HTB_(b, a) = complex<double>(1.0 * Parameters_.t_hopping, 0.0) * phasex;
                HTB_(a, b) = conj(HTB_(b, a));
            }
        }

        // * +y direction Neighbor
        if (Coordinates_.indy(l) == (Coordinates_.ly_ - 1))
        {
            phasex = one_complex;
            phasey = exp(iota_complex * 2.0 * (1.0 * my) * PI / (1.0 * Parameters_.TBC_cellsY));
        }
        else
        {
            phasex = one_complex;
            phasey = one_complex;
        }
        m = Coordinates_.neigh(l, 2);
        for (int spin = 0; spin < 2; spin++)
        {

            a = l + ns_ * spin;
            b = m + ns_ * spin;
            assert(a != b);
            if (a != b)
            {

                HTB_(b, a) = complex<double>(1.0 * Parameters_.t_hopping, 0.0) * phasey;

                HTB_(a, b) = conj(HTB_(b, a));
            }
        }
    }

} // ----------

void Hamiltonian::HTBClusterCreate()
{

    int ns = (Parameters_.lx_cluster) * (Parameters_.ly_cluster);

    complex<double> phasex, phasey;
    int l, m, a, b;

    HTBCluster_.fill(0.0);

    for (l = 0; l < ns; l++)
    {

        // * +x direction Neighbor
        if (CoordinatesCluster_.indx(l) == (CoordinatesCluster_.lx_ - 1))
        {
            phasex = one_complex;
            phasey = one_complex;
        }
        else
        {
            phasex = one_complex;
            phasey = one_complex;
        }
        m = CoordinatesCluster_.neigh(l, 0);
        for (int spin = 0; spin < 2; spin++)
        {

            a = l + ns * spin;
            b = m + ns * spin;
            assert(a != b);
            if (a != b)
            {

                HTBCluster_(a, b) = complex<double>(1.0 * Parameters_.t_hopping, 0.0) * phasex;
                HTBCluster_(b, a) = conj(HTBCluster_(a, b));
            }
        }

        // * +y direction Neighbor
        if (CoordinatesCluster_.indy(l) == (CoordinatesCluster_.ly_ - 1))
        {
            phasex = one_complex;
            phasey = one_complex;
        }
        else
        {
            phasex = one_complex;
            phasey = one_complex;
        }
        m = CoordinatesCluster_.neigh(l, 2);
        for (int spin = 0; spin < 2; spin++)
        {

            a = l + ns * spin;
            b = m + ns * spin;
            assert(a != b);
            if (a != b)
            {

                HTBCluster_(a, b) = complex<double>(1.0 * Parameters_.t_hopping, 0.0) * phasey;
                HTBCluster_(b, a) = conj(HTBCluster_(a, b));
            }
        }
    }
    // HTBCluster_.print();

} // ----------

void Hamiltonian::Hoppings()
{
    //DOES SOMETHING EXACT i.e NOTHING :)

} // ----------

void Hamiltonian::copy_eigs(int i)
{

    int space = 2 * ns_;

    if (i == 0)
    {
        for (int j = 0; j < space; j++)
        {
            eigs_[j] = eigs_saved_[j];
        }
    }
    else
    {
        for (int j = 0; j < space; j++)
        {
            eigs_saved_[j] = eigs_[j];
        }
    }
}

void Hamiltonian::copy_eigs_Cluster(int i)
{

    int ns = (Parameters_.lx_cluster) * (Parameters_.ly_cluster);
    int space = 2 * ns;

    if (i == 0)
    {
        for (int j = 0; j < space; j++)
        {
            eigsCluster_[j] = eigsCluster_saved_[j];
        }
    }
    else
    {
        for (int j = 0; j < space; j++)
        {
            eigsCluster_saved_[j] = eigsCluster_[j];
        }
    }
}

#endif
