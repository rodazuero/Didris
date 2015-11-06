// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/random.hpp>
#include <boost/random/variate_generator.hpp>
#include <ctime>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
using namespace arma;
using namespace Rcpp;

//I. Solving the optimal level of a_{t+1}. First I define a c-structure of
//parameters for the solver


//===================================================================================//
//Defining the structure of parameters for the asset solver
//===================================================================================//
struct asset_my_params {double Astruct_b;  double Astruct_wl; double Astruct_wh;  double Astruct_Ph; double Astruct_h; double Astruct_eeta;  double Astruct_i; double Astruct_ppsi; double Astruct_ttheta;  double Astruct_bbeta;};


//Defining the function to which we will define the roots
double F_root_a (double a, void *params){
    //0. Assigning the parameters
    struct asset_my_params *p=(struct asset_my_params *) params;
    double b=p->Astruct_b;
    double wl=p->Astruct_wl;
    double wh=p->Astruct_wh;
    double Ph=p->Astruct_Ph;
    double h=p->Astruct_h;
    double eeta=p->Astruct_eeta ;
    double i=p->Astruct_i ;
    double ppsi=p->Astruct_ppsi;
    double ttheta=p->Astruct_ttheta;
    double bbeta=p->Astruct_bbeta;
    //1. Getting the function residual to obtain the root
    double residual=pow(b-a+wl*(1-h)-Ph*h,-ppsi)-
    bbeta*(1+i)*(pow(h*ttheta,eeta)*pow((wh+(1+i)*a),-ppsi)+
    (1-pow(h*ttheta,eeta))*pow(wl+(1+i)*a,-ppsi));
    return(residual);
}

//=========================
//Solver for the asset
//=========================
// [[Rcpp::export]]
long double Assetsolver(double b, double wl, double wh, double Ph, double h,
                        double eeta,double i, double ppsi, double ttheta,
                        double bbeta, double Amax){
    int status;
    int iter=0, max_iter=100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r=0, r_expected=22;
    double x_lo=-wl/(1+i)+1.0e-3, x_hi=b+wl*(1-h)-Ph*h-1.0e-3;
    gsl_function F; struct asset_my_params params={b, wl, wh, Ph, h, eeta, i,
        ppsi, ttheta, bbeta};
    F.function= &F_root_a;
    F.params= &params;
    cout << F.params << endl;
    T = gsl_root_fsolver_brent;
    s=gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    //printf ("using %s method /n", gsl_root_fsolver_name (s));
    //printf ("%5s [%9s, %9s] %9s %10s %9s/n",
    //   "iter" , "lower" , "upper" , "root" , "err" , "err(est)");
    //cout << " whe are here 6" << endl;
    do{
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi,0, 0.001);
        
        //if (status == GSL_SUCCESS)
         //printf ("Converged:\n");
        //printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
        //iter, x_lo, x_hi,r, r - r_expected, x_hi - x_lo);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);
    //Borrowing constraint
    if (r<=Amax){
        r=Amax;
    }
    return(r);
}
//=========================
//Utility function 
//=========================
// [[Rcpp::export]]
long double Hopt(vec par, double b, double ttheta,double wl, double wh, double i,
                 double Ph, double Amax){
    //0. Initializing utility if education or not
    double U1=0;
    double U0=0;

    
    //1. Loading parameters
    double bbeta=par[0];
    double eeta=par[1];
    double ppsi=par[2];
    
    //2. Get optimal assets holding in both scenarios: educ or not
    double a0=Assetsolver(b,wl, wh,  Ph,  0, eeta, i,  ppsi,  ttheta,
                          bbeta,  Amax);
    double a1=Assetsolver(b,wl, wh,  Ph,  1, eeta, i,  ppsi,  ttheta,
                          bbeta,  Amax);
    //3. Comparing utilities in both cases
    U0=pow(b-a0+wl,1-ppsi)/(1-ppsi)+bbeta*pow(wl+(1+i)*a0,1-ppsi)/(1-ppsi);
    U1=pow(b-a1-Ph,1-ppsi)/(1-ppsi)+bbeta*(
        pow(ttheta,eeta)*pow(wh+(1+i)*a1,1-ppsi)/(1-ppsi)+
        (1-pow(ttheta,eeta))*pow(wl+(1+i)*a1,1-ppsi)/(1-ppsi));
    return(U1>=U0);
}


// [[Rcpp::export]]
double F_root_ab (double a, double wl, double h, double Ph,
                  double ppsi, double bbeta, double ttheta,
                  double wh, double eeta, double b, double i){
    //0. Assigning the parameters

    //1. Getting the function residual to obtain the root
    double residual=pow(b-a+wl*(1-h)-Ph*h,-ppsi)-
    bbeta*(1+i)*(pow(h*ttheta,eeta)*pow((wh+(1+i)*a),-ppsi)+
                 (1-pow(h*ttheta,eeta))*pow(wl+(1+i)*a,-ppsi));
    return(residual);
}

