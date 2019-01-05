#include <cmath>
#include <vector>

using namespace std;

// ==============================================================================
/**
Class handling plasma parameters e.g. temperature, density, composition, magnetic fields etc.
Uses SI units; temperatures in Jouls.
*/
// ==============================================================================

class plasma : public
{

public:	

  // Plasma composition
  vector<double> mass;
  vector<double> charge;
  int electronIndex;

  void temperature(double* val, const double* x, const int indexSpecies) const;
  void density(double* val, const double* x, const int indexSpecies) const;

private:
  double _electron;
  double _massElectron;
  double _atomicMassUnit;

  boolean _rho_tor_present;
  int     _rho_tor_index;

  boolean _velocity_present;
  int     _velocity_index;

  boolean _energy_present;
  int     _energy_index;
  char*   _energy_unit;

  double _temperature0;
  vector<double> _rho_tor_grid1d;
  vector<double> _temperature_grid1d;

  double _densityElectrons0;
  vector<double> _density0;

  double plasma::interp1d(vector<double> x, vector<double> y, double xout) const;
};

// End of Class

// ==============================================================================
// IMPLEMENTATIONS
// ==============================================================================
//
// The constructor
plasma::plasma(xmllib XMLparams){

  char* varnames=getXMLvalue(XMLparams,'variableNames');
  for (int j=0; j < size(varnames)   ; j++ ) {
    if (varnames[j]=='rho_tor'){
      _rho_tor_present = true;
      _rho_tor_index   = j;
    }
    if (varnames[j]=='velocity'){
      _velocity_present = true;
      _velocity_index   = j;
    }
    if (varnames[j]=='energy'){
      _energy_present = true;
      _energy_index   = j;
      //_energy_unit    = 'joul'
    }
  }

  _electron = 1.60217656535e-19;
  _atomicMassUnit = 1.66053892173e-27;
  _massElectron = 9.1093829140e-31;

  /////////////////////////////////////////////
  // BEGIN: Default values
  /////////////////////////////////////////////
  _temperature0 = 20e3 * _electron;

  mass.clear();
  charge.clear();
  // Electron
  charge.push_back( -1.0 * _electron );
  mass.push_back(    1.0 * _massElectron );
  electronIndex = _mass.size;

  // Deuterium
  mass.push_back(    2.01410178 * _atomicMassUnit );
  charge.push_back(  1.0        * _electron );
  // Tritium
  mass.push_back(    3.0160492  * _atomicMassUnit );
  charge.push_back(  1.0        * _electron );
  // He3
  mass.push_back(    3.0160293  * _atomicMassUnit );
  charge.push_back(  2.0        * _electron );
  // He4
  mass.push_back(    4.0026022  * _atomicMassUnit );
  charge.push_back(  2.0        * _electron );
  // Carbon12
  mass.push_back(   12.01078    * _atomicMassUnit );
  charge.push_back(  6.0        * _electron );

  double ne0 = 1e20;
  density0.clear();
  density0.push_back( ne0              );  // Electrons
  density0.push_back( ne0 * 0.40 / 1.0 );  // Deuterium
  density0.push_back( ne0 * 0.40 / 1.0 );  // Tritium
  density0.push_back( ne0 * 0.03 / 3.0 );  // He3
  density0.push_back( ne0 * 0.12 / 4.0 );  // He4
  density0.push_back( ne0 * 0.05 / 6.0 );  // Carbon

  _rho_tor_grid1d.clear();
  _rho_tor_grid1d.push_back(0.0);
  _rho_tor_grid1d.push_back(1.0);

  _temperature_grid1d.clear();
  _temperature_grid1d.push_back( _temperature0 );
  _temperature_grid1d.push_back( _temperature0 );

  /////////////////////////////////////////////
  // END: Default values
  /////////////////////////////////////////////
}


void plasma::temperature(double* val, const double* x, const int indexSpecies) const{
  if (_rho_tor_present){
    val[0] = interp1d(_rho_tor_grid1d , _temperature_grid1d , x[_rho_tor_index]);
  } else {
    val[0] = _temperature0;
  }
}

void plasma::density(double* val, const double* x, const int indexSpecies) const{
  val[0] = _density0[indexSpecies];
}

double plasma::interp1d(vector<double> x, vector<double> y, double xout) const {
  return y[0];
}
