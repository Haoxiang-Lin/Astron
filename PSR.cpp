#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

//***************This is the space allowed to modify***************
//Pulsar Parameters
const double        yr2sec      = 31536000.;
const double        T_SNR       = 957. * yr2sec;
const double        n_ob        = 2.51;
const double        T_Ch_ob     = 1240. * yr2sec;
const double        m_ob        = 10.23;
const double        period      = 0.033;
//allowed range of n, T_Ch and m
const double        n_ob_min    = n_ob * 0.8;
const double        n_ob_max    = n_ob * 1.2;
const double        T_Ch_ob_min = T_Ch_ob * 0.8;
const double        T_Ch_ob_max = T_Ch_ob * 1.2;
const double        m_ob_min    = m_ob * 0.8;
const double        m_ob_max    = m_ob * 1.2;
//range of T_SNR and Gamma
const double        Tmin        = 0.1;
const double        Tmax        = 2.0;
const double        rmin        = -1.0;
const double        rmax        = 2.0;
//step of T_SNR and Gamma
const int           Tstep       = 190;
const int           rstep       = 300;
//range of B_* and dM(t0)
const double        Bmin        = 10.;
const double        Bmax        = 13.;
const double        dMmin       = 23.;
const double        dMmax       = 30.;
//step of of B_* and dM(t0)
const int           Bstep       = 20;
const int           dMstep      = 40;
//*****************************************************************

//********These are the constants used in this calculation*********
const long double   pi          = 3.14159265358979323846264338327950;
const float         G           = 6.673e-8F;
const double        c           = 2.99792458e010;
const float         M_Sun       = 1.9891e33F;

const double        M0          = 1.4 * M_Sun;
const double        R0          = pow(10., 6);
const double        I           = 0.4 * M0 * R0 * R0;
const double        t0          = 300.;
const double        w           = 2 * pi / period;
//*****************************************************************

int main()
{
    //Output judgement on a specific point 
    /*
    bool PSR_BOOL(double r, double T_SNR);
    cout << PSR_BOOL(-1, 957) << endl;
    */

    void PSR_CALL(string file_name);

    string file_name;
    cout << "Enter the name of the file to save the data:      ";
    cin >> file_name;

    PSR_CALL(file_name);

    system("pause");
    return 0;
}

//function to print results
void PSR_CALL(string file_name)
{
    //set up output stream
    char xpause;
    ofstream outdata;
    outdata.open(file_name.c_str(), ios::out);
    if (outdata.fail()){
        cout << " Unable to open file --- Terminating calculation" 
            << "\n\n"
            << "Enter any character to quit: ";
        cin  >> xpause;   
        exit(1);
    }

    //declare functions
    void linspace(double *begin, double *end, const double min, const double max, const int step);
    bool PSR_BOOL(double r, double T_SNR);

    //array of T_SNR
    double t_[Tstep + 1] = {0};
    double *pt = t_;
    linspace(pt, pt + Tstep + 1, Tmin * T_SNR, Tmax * T_SNR, Tstep);

    //array of Gamma
    double r_[rstep + 1] = {0};
    double *pr = r_;
    linspace(pr, pr + rstep + 1, rmin, rmax, rstep);

    //print output on the prompt and to the file: (r, T/T_SNR, bool)
    cout    << fixed << showpoint << right << setprecision(6);
    outdata << fixed << showpoint << right << setprecision(6);
    for (int i = 0; i != rstep + 1; ++i)
        for (int j = 0; j != Tstep + 1; ++j) {
            cout    << setw(15) << r_[i] << setw(15) << t_[j]/T_SNR << setw(15) << PSR_BOOL(r_[i], t_[j]) << endl;
            outdata << setw(15) << r_[i] << setw(15) << t_[j]/T_SNR << setw(15) << PSR_BOOL(r_[i], t_[j]) << endl;
        }
}

//Function to judge if (r, T_SNR) satisfies condition
bool PSR_BOOL(double r, double T_SNR)
{
    //declare functions
    void PSR_CALCULATE(double r, double T_SNR, double B, double dM_t0, double& n, double& T_Ch, double& m);
    void logspace(double *begin, double *end, const double min, const double max, const int step);

    //define n, T_Ch and m
    double n, T_Ch, m;

    //array of B_*
    double B_[Bstep + 1]   = {0};
    double *pB  = B_;
    logspace(pB,  pB + Bstep + 1,   Bmin,   Bmax,   Bstep);

    //array of dM(t0)
    double dM_[dMstep + 1] = {0};
    double *pdM = dM_;
    logspace(dM_, dM_ + dMstep + 1, dMmin,  dMmax,  dMstep);

    //calculate and make judgement
    for (int i = 0; i != Bstep + 1; ++i)
        for (int j = 0; j != dMstep + 1; ++j) {
            PSR_CALCULATE(r, T_SNR, B_[i], dM_[j], n, T_Ch, m);
            if ((n <= n_ob_max) && (n >= n_ob_min) && (T_Ch <= T_Ch_ob_max) && (T_Ch >= T_Ch_ob_min) 
                && (m <= m_ob_max) && (m >= m_ob_min))
                return true;
        }
    return false;
}

//function to calculate n and T_Ch, with definition of equations
void PSR_CALCULATE(double r, double T_SNR, double B, double dM_t0, double& n, double& T_Ch, double& m)
{
    //magnetic moment
    double mu    = B * pow(R0, 3);
    //mass accretion rate
    double dM_t  = dM_t0 * pow((T_SNR/t0), -1.25);
    //magnetosphere radius
    double R_m   = pow( pow(mu, 2)/(dM_t * sqrt(2 * G * M0) ), 2/7.);
    //light cylindrical radius
    double R_LC  = c / w;
    //Keplerian angular velocity
    double w_K   = sqrt(G * M0 / pow(R_m, 3));
    //magnetic spin-down torque
    double T_mag = - 0.5 * pow(mu, 2) * pow(R_m, -2) * pow(R_LC, -1);
    //propeller spin-down torque
    double T_p   = - dM_t * pow(R_m, 2) * w_K * pow(w/w_K, r);
    //first derivative of Omega
    double dw    = (T_mag + T_p)/I;
    //second derivative of Omega
    double ddw   = (- 5/(7*T_SNR) + dw/w) * T_mag/I + ((-30+15*r)/(28*T_SNR) + r*dw/w) * T_p/I;
    //third derivative of Omega
    double dddw  = pow((- 5/(7*T_SNR) + dw/w), 2) * T_mag/I + pow(((-30+15*r)/(28*T_SNR) + r*dw/w), 2) * T_p/I 
        + (5/(7*pow(T_SNR, 2)) + (ddw*w - pow(dw, 2))/pow(w, 2)) * T_mag/I
        + ((30-15*r)/(28*pow(T_SNR, 2)) + r*(ddw*w - pow(dw, 2))/pow(w, 2)) * T_p/I;
    //braking index
    n    = ddw*w/pow(dw, 2);
    //characteristic age
    T_Ch = -w/(2*dw);
    //secondary braking index
    m    = dddw*pow(w, 2)/pow(dw, 3);
}

//function to generate log10 array
void logspace(double *begin, double *end, const double min, const double max, const int step)
{
    double* previous = begin;
    *begin = pow(10., min);
    begin ++;
    while (begin != end) {
        *begin = *previous * pow(10., (max - min) / double(step));
        previous++;
        begin++;
    } 
}

//function to generate linear array
void linspace(double *begin, double *end, const double min, const double max, const int step)
{
    double* previous = begin;
    *begin = min;
    begin ++;
    while (begin != end) {
        *begin = *previous + (max - min) / double(step);
        previous++;
        begin++;
    } 
}

