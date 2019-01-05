#include <cmath>
#include <vector>
#include <plasma.hpp>

using namespace std;

// ==============================================================================
/**
Class handling coloumb collisions against a Maxwellian background plasma.
Includes the Chandrasekar/Stix coefficients alpha, beta and gamma.
Uses SI units; temperatures in Jouls.
*/
// ==============================================================================

class coloumbCollisionsOnMaxwellian : public
{

public:	

  double coloumbLogarithm(const double* x, const int indexBackgroundSpecies) const;

  void alpha(double* val, const double* x, const int indexSpecies) const;
  void alphaOnSingleSpecies(double* val, const double* x, const int indexSpecies, const int indexBackgroundSpecies) const;

  void beta(double* val, const double* x, const int indexSpecies) const;
  void betaOnSingleSpecies(double* val, const double* x, const int indexSpecies, const int indexBackgroundSpecies) const;

  void gamma(double* val, const double* x, const int indexSpecies) const;
  void gammaOnSingleSpecies(double* val, const double* x, const int indexSpecies, const int indexBackgroundSpecies) const;

  void AlphaBetaGamma_OnSingleSpecies(double* alpha, double* beta, double* gamma, const double* vel, const double* x, const int indexSpecies, const int indexBackgroundSpecies) const;
  void AlphaBetaGamma_OnSingleSpecies(double* alpha, double* beta, double* gamma, const double* x, const int indexSpecies, const int indexBackgroundSpecies) const;

private:
  // Plasma composition
  plasma* _plasma;
  double velocity(double* x, int indexSpecies) const;
  double Cf_Stix(double* x, int indexSpecies, int indexBackgroundSpecies) const;
  double lf_Stix(double* x, int indexSpecies, int indexBackgroundSpecies) const;
  double G_Stix(double X) const;
};

// ==============================================================================
//  CONSTRUCTOR
// ==============================================================================

coloumbCollisionsOnMaxwellian::coloumbCollisionsOnMaxwellian(plasma* plasmaIn, int indexSpecies){
  _plasma = plasmaIn;
}

// ==============================================================================
//  PUBLIC
// ==============================================================================

double coloumbCollisionsOnMaxwellian::coloumbLogarithm(const double* x, const int indexBackgroundSpecies) const {
  double Temp;
  double Dens;
  _plasma.temperature(&Temp, x, indexBackgroundSpecies);
  _plasma.electronDensity(&Dens, x, indexBackgroundSpecies);
  if (indexBackgroundSpecies==_plasma.electronIndex){
    return 15.2 - 0.5*log(Dens/1e20) + log(Temp/_electron);
  } else {
    return 17.3 - 0.5*log(Dens/1e20) + 1.5*log(Temp/_electron);
  }
}

// ==============================================================================
//  PRIVATE
// ==============================================================================

double coloumbCollisionsOnMaxwellian::velocity(double* x, int indexSpecies) const{
  if (plasma._velocity_present){
    return x[plasma._velocity_index];
  }
  if (plasma._energy_present){
    return sqrt( 2.0 * x[plasma._energy_index] / _plasma.mass[indexSpecies] );
  }
  return 0.0;
}

double coloumbCollisionsOnMaxwellian::lf_Stix(double* x, int indexSpecies, int indexBackgroundSpecies) const{
  double Temp;
  _plasma.temperature(&Temp , x , indexBackgroundSpecies);
  return sqrt(mass[indexSpecies] / (2.0*Temp) );
}

double coloumbCollisionsOnMaxwellian::Cf_Stix(double* x, int indexSpecies, int indexBackgroundSpecies) const{
  return (0.5/(M_PI*pow(8.85418782e-12,2))) *  pow(charge[indexSpecies]*charge[indexBackgroundSpecies]/mass[indexSpecies],2) * coloumbCollisionsOnMaxwellian(x,indexBackgroundSpecies);
}

double coloumbCollisionsOnMaxwellian::G_Stix(double X) const{
  return (erfc(X) + exp(-X*X)*2.0/sqrt(M_PI)) / (2.0*X);
}


double coloumbCollisionsOnMaxwellian::alphaOnSingleSpecies(const double* x, const int indexSpecies, const int indexBackgroundSpecies) const{
  double vel = velocity( x , indexSpecies );
  double lf  = lf_Stix( x , indexSpecies, indexBackgroundSpecies);
  double G   = G_Stix( vel * lf );
  double Cf  = Cf_Stix(j,jf);
  return lf*lf*Cf * (1.0-mass[indexSpecies]/mass[indexBackgroundSpecies]) * G;
}

double coloumbCollisionsOnMaxwellian::betaOnSingleSpecies(const double* x, const int indexSpecies, const int indexBackgroundSpecies) const{
  double vel = velocity( x , indexSpecies );
  double lf  = lf_Stix( x , indexSpecies, indexBackgroundSpecies);
  double G   = G_Stix( vel * lf );
  double Cf  = Cf_Stix(j,jf);
  return (Cf / vel) * G;
}

double coloumbCollisionsOnMaxwellian::gammaOnSingleSpecies(const double* x, const int indexSpecies, const int indexBackgroundSpecies) const{
  double vel = velocity( x , indexSpecies );
  double lf  = lf_Stix( x , indexSpecies, indexBackgroundSpecies);
  double G   = G_Stix( vel * lf );
  double phi = erfc(x);
  double Cf  = Cf_Stix(j,jf);
  return (Cf / vel) * (phi-G);
}

void coloumbCollisionsOnMaxwellian::AlphaBetaGamma_OnSingleSpecies(double* alpha, double* beta, double* gamma, const double* vel, const double* x, const int indexSpecies, const int indexBackgroundSpecies) const{
  double vel = velocity( x , indexSpecies );
  double lf  = lf_Stix( x , indexSpecies, indexBackgroundSpecies);
  double G   = G_Stix( vel * lf );
  double phi = erfc(x);
  double Cf  = Cf_Stix(j,jf);
  alpha[0] = Cf * pow(lf,2) * (1.0-mass[indexSpecies]/mass[indexBackgroundSpecies]) * G;
  beta[0]  = (Cf / vel) * (phi-G);
  gamma[0] = (Cf / vel) * (phi-G);
}

void coloumbCollisionsOnMaxwellian::AlphaBetaGamma_OnSingleSpecies(double* alpha, double* beta, double* gamma, const double* x, const int indexSpecies, const int indexBackgroundSpecies) const{
  double vel = velocity( x , indexSpecies );
  AlphaBetaGamma_OnSingleSpecies(alpha, beta, gamma, vel, x, indexSpecies, j );
}

double coloumbCollisionsOnMaxwellian::alpha(const double* x, const int indexSpecies) const{
  double v=0;
  for (int j; j<size(_plasma.mass); j++){
    v += alphaOnSingleSpecies(&vNew , x , indexSpecies , j);
  }
  return v;
}

double coloumbCollisionsOnMaxwellian::beta(const double* x, const int indexSpecies) const{
  double v=0;
  for (int j; j<size(_plasma.mass); j++){
    v += betaOnSingleSpecies(&vNew , x , indexSpecies , j);
  }
  return v;
}

double coloumbCollisionsOnMaxwellian::gamma(const double* x, const int indexSpecies) const{
  double v=0;
  for (int j; j<size(_plasma.mass); j++){
    v += gammaOnSingleSpecies(&vNew , x , indexSpecies , j);
  }
  return v;
}

void coloumbCollisionsOnMaxwellian::AlphaBetaGamma(double* alpha, double* beta, double* gamma, const double* x, const int indexSpecies) const{
  double alphaNew;
  double betaNew;
  double gammaNew;
  alpha[0]=0.0;
  beta[0]=0.0;
  gamma[0]=0.0;
  double vel = velocity( x , indexSpecies );
  for (int jBackground; jBackground<size(_plasma.mass); jBackground++){
    AlphaBetaGamma_OnSingleSpecies(alphaNew, betaNew, gammaNew, vel, x, indexSpecies, jBackground );
    alpha[0] += alphaNew;
    beta[0]  += betaNew;
    gamma[0] += gammaNew;
  }
}
