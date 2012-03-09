#include "Delay.h"

// Static members
const unsigned Delay::n_par = 2 ;

// Element information
ItemInfo Delay::einfo =
{
  "delay",
  "Ideal bi-directional circuit delay",
  "Chris Saunders",
  DEFAULT_ADDRESS"elements/Delay.h.html"
};

// Parameter information
ParmInfo Delay::pinfo[] =
{
  {"delay", "Time Delay",TR_DOUBLE, true},
  {"res" , "Internal conditioning resistance", TR_DOUBLE, true}
};


Delay::Delay(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  // Value required
  paramvalue[0] = &(delay=0) ;
  paramvalue[1] = &(res=50);

  // Set flags
  setFlags(NONLINEAR | MULTI_REF | TR_TIME_DOMAIN);
}

void Delay::getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0)); // Input port terminal
  term_list.push_back(getTerminal(1)); // Input port local reference terminal
  local_ref_vec.push_back(1); // Local reference index

  term_list.push_back(getTerminal(2)); // Output port terminal
  term_list.push_back(getTerminal(3)); // Output port local reference terminal
  local_ref_vec.push_back(3); // Local reference index
}

void Delay::init() throw(string&)
{
  if(!isSet(&delay) && (delay<0.0) )
  {
    throw(" Delay must be positive! ") ;
  }
  else
  {
    // Set the number of terminals
    setNumTerms(4) ;

    // Set number of states
    setNumberOfStates(2);

    DenseIntVector var(2,0); // Create 4 state variables for I1, I2, g1, g2
    var[0]=0;
    var[1]=1;
    DenseIntVector novar;
    DenseIntVector delay_vec(2,0);  //Create delayed state variables for g1, g2
    delay_vec[0]=0;
    delay_vec[1]=1;
    DenseDoubleVector delay_t_vec(2,0);  //Set both time delays to the declared value
    delay_t_vec[0]= delay;
    delay_t_vec[1]= delay;
    initializeAD(var, novar, novar, delay_vec, delay_t_vec) ;
    temp_eval_counter = 0;
  }
}

void Delay::eval(AD * x, AD * effort, AD * flow)
{
  // State Variables:
  // effort[0] - Port 1 voltage
  // effort[1] - Port 2 voltage
  // flow[0] - Port 1 input current
  // flow[1] - Port 2 input current
  // x[0] - g1
  // x[1] - g2
  // x[2] - g1(t - delay)
  // x[3] - g2(t - delay)

  /*
  AD temp_res = 0.5 * res * ((x[3] +x [2]) / (x[3] - x[2] + 1e-15));
  cout << temp_res << endl;
  */

  //cout << this->getCurrentTime() <<endl;

  flow[0] = 0.5 * (x[1] - x[2]);
  flow[1] = 0.5 * (x[0] - x[3]);

  effort[0] = 0.5 * (x[1] + x[2]) * res;
  effort[1] = 0.5 * (x[0] + x[3]) * res;

  /*
  // Coeffs to model attenuation
  flow[0] = 0.5 * (x[1] - 0.8 * x[2]);
  flow[1] = 0.5 * (x[0] - 0.8 * x[3]);

  effort[0] = 0.5 * (x[1] + 0.8 * x[2]) * res;
  effort[1] = 0.5 * (x[0] + 0.8 * x[3]) * res;
  */

  //  temp_eval_counter++;
  //  cout << "Eval routine iterations: " << temp_eval_counter << endl;

  /*
  AD temp = (flow[0] - flow[1])/(x[0] + x[1]);
  cout << temp << endl;
  */

  //temp_res = 0.5 * (effort[0] + effort[1])/(flow[0] + flow[1]);

  /*
  flow[0] = x[0];
  flow[1] = x[1];
  */

  /*
  effort[0] = (x[4] + x[0])*res;
  effort[1] = (x[5] + x[1])*res;

  x[2] = (effort[1]/res) + x[1];
  x[3] = (effort[0]/res) + x[0];
  */

  /*
  effort[0] = (x[4] + x[0])*res + (x[3] - x[2] - 2*x[0])*10;
  effort[1] = (x[5] + x[1])*res + (x[2] - x[3] - 2*x[1])*10;
  */

  /*
  condassign(x[2],1,(effort[1]/res) + flow[1],(effort[1]/res) + flow[1]);
  condassign(x[3],1,(effort[0]/res) + flow[0],(effort[0]/res) + flow[0]);
  */

  //cout << "v1= " << effort[0] << " v2= " << effort[1] << " i1= " << flow[0] << " i2= " << flow[1] << " g1= " << x[0] << " g2= " << x[1] << " g1d= " << x[2] << " g2d= " << x[3] << endl;
}
