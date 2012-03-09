// Contains routines for analysis of IBIS model based on
// 'A Novel Extraction Method of Analog
// SPICE Behavioral Model from IBIS  Model'
// Hwan-Mok Jung, Chang-Gene Woo, Pyung Choi,
// Jong-Hwa Kwon, Jae-Hoon Yun
//
//			o Ramp Voltage with slope 1
//    			|
//              |-------------|
//              |             |
// Input  o-----|             |----o Buffer output
//              |             |
//              |--------------
//                     |
//                     |
//                     o
//                 Reference
//
//
//   Authors:
//
//   Nikhil Mahajan,Ambrish Varma and Rubina Ahmed
//

#ifndef Ibis_h
#define Ibis_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class Ibis : public ADInterface
{
	double tCrossed; //globals to determine whether Vin has crossed
	//Vth
	adouble   t_call;

	public:

	Ibis(const string& iname);
	~Ibis() {}

	static const char* getNetlistName()
	{
    return einfo.name;
	}

	// Do some local initialization
	virtual void init() throw(string&);

	private:

	virtual void eval(AD * x, AD * effort, AD * flow);

	// Element information
	static ItemInfo einfo;

	// Number of parameters of this element
	static const unsigned n_par;

	// Parameter variables
	string Ibis_file;
	double Vcc;

	// Parameter information
	static ParmInfo pinfo[];

	static int n_eval_calls,r_f_old,toggleME;
	static adouble x0_old;
};

#endif
