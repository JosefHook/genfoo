//
// Define system parameters
//
#pragma once
#ifndef RF1DPARAMS_H
#define RF1DPARAMS_H


#include <cmath>

//////////////////////////////////////////
// SI Static parameters
//////////////////////////////////////////
long double Ze = -1.0;               // Electron number
long double Zp = 1.0;                // proton number
long double Z = 1.0;                 // testparticle number
long double np = 1e20;               // Density protons             [1/L^3]     
long double ne = np;                 // Density electrons           [1/L^3]     
long double nh = 0.02*np;            // Density of heated particles 
long double kb = 1.3807e-23;         // Boltzman [J/K] 
long double Te = 20e3*1.1604505e4;   // Electron temperature       [eV]*[K/eV] = [K]
long double Tp = 20e3*1.1604505e4;   // Proton temperature         [eV]*[K/eV] = [K]
long double T =  20e3*1.1604505e4;   // test particle temperature  [eV]*[K/eV] = [K] 
long double mp = 1.6726e-27;         // Mass proton
long double mee = 9.1094e-31;      // Mass electron



long double m = mp*2.0;                // Mass testparticle
long double q = 1.6022e-19;          // elemeentary charge
long double eps0 = 8.8542e-12;       // permativity
long double Power = 0.2e6;           // Typical RF power for JET is 32MW
long double K=Power/(3.0*m*nh);

long double CONV=1.173702880724304e-015; //Conversion factor, from CGS to SI units
long double loglambda = std::log( CONV*3/(2*Z*Zp*std::pow(q,3))*sqrt( std::pow(kb*T,3)/(M_PI*np)) );
long double ts = 3*std::pow(2.0*M_PI,3.0/2.0)*std::pow(eps0,2)*m*std::pow(kb*Te,3.0/2.0)/( std::pow(Z,2)*std::pow(q,4)*sqrt(mee)*ne*loglambda  );
  
long double lp = sqrt(mp/(2.0*kb*Tp));


long double le = sqrt(    2.0*(mee/kb)/Te  );//sqrt(mee/(2*kb*Te));




long double vth = 1/lp;


long double lelp = le/lp; //sqrt( (mee*Tp)/(mp*Te) );


//
// Size of parameeters

// Plasma temperature: 2.3209e+08
// Plasma density    : 2e+18
// RF power          : 200000
// K           : 9.96453e+12
// K/vth           : 5.09049e+06
// Cp                : 2.65404e+20
// Ce                : 2.65404e+20
// np                : 1e+20
// std::pow(Zp,2)    : 1
// std::pow(Z,2)     : 1
// std::pow(q,4)     : 6.58972e-76
// log lambda        : 22.2005
// std::pow(4*pi*eps0*m,2) : 1.38536e-73
// eps0 : 8.8542e-12
// m : 3.3452e-27
// vth : 1.95748e+06
// p1                : 4.78955e+07
// p2                : 22.2005
// p3                : 1.80954e-08
// p4                : 1.80954e+12
// p5                : 1.80954e+12


long double p1 = q/m;
long double p2 = std::pow(Z,2)*loglambda;
long double p3 = q/eps0;
long double p4 = np*q/eps0;
long double p5 = ne*q/eps0;

// Old 
//long double Cp  = 8*M_PI*np*std::pow(Zp,2)*std::pow(Z,2)*std::pow(q,4)*loglambda/( std::pow(4*M_PI*eps0*m,2));

long double Cp = (8.0*M_PI*p2*std::pow(Zp,2))*(p1*p4)*(p1*p3)/( std::pow(4.0*M_PI,2));

// Oold
//long double Ce  = 8*M_PI*ne*std::pow(Ze,2)*std::pow(Z,2)*std::pow(q,4)*loglambda/( std::pow(4*M_PI*eps0*m,2));

long double Ce = (8.0*M_PI*p2*std::pow(Ze,2))*(p1*p5)*(p1*p3)/( std::pow(4.0*M_PI,2));

#endif
