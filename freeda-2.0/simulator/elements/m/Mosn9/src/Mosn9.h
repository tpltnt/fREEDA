// Philips MOS9 model
//                Drain 2
//                  o
//                  |
//                  |
//              |---+
//              |
// Gate 1 o-----|------o 4 Bulk
//              |
//              |---+
//                  |
//                  |
//                  o
//               Source 3
//
//	  Author: Kuldip

#ifndef Mosn9_h
#define Mosn9_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class Mosn9:public ADInterface
{
	public:

	Mosn9(const string& iname);
	~Mosn9() {}

	//Element informations
	static ItemInfo einfo;

	static const char* getNetlistName()
	{
		return einfo.name;
	}

	// do some local initialization
	virtual void init() throw(string&);

	private:

	//implement the eval function
	virtual void eval(AD * x, AD * effort, AD * flow);

	//number of parameters of this element
	static const unsigned n_par;

	// parameter variables
	double LEVEL,VT0,K0,K,PHIB,VSBX,BET,THE1,THE2,THE3,GAM1;
	double ETADS,ALP,VP,GAM00,ETAGAM,M0,ETAM,PHIT,ZET1,VSBT,A1,A2,A3;
	double COX,CGDO,CGSO,NT ,NFMOD,NF,NFA,NFB,NFC,TOX,MULT;
	double T0,BOLTZMAN,Q,EPIOX;

  // Parameter variables
  //Parameters of the geometrical model
  double LER,WER,LVAR,LAP,WVAR,WOT,TR,VT0R,STVT0;
  double SLVT0,SL2VT0,SL3VT0,SWVT0,K0R,SLK0,SL2K0,SWK0,KR,SLK;
  double SL2K,SWK,PHIBR,VSBXR,SLVSBX,SWVSBX,BETSQ,ETABET,LP1,FBET1;
  double LP2,FBET2,THE1R,STTHE1R,SLTHE1R,STLTHE1,GTHE1,SWTHE1,WDOG,FTHE1;
  double THE2R,STTHE2R,SLTHE2R,STLTHE2,SWTHE2,THE3R,STTHE3R,SLTHE3R,STLTHE3,SWTHE3;
  double GAM1R,SLGAM1,SWGAM1,ETADSR,ALPR,ETAALP,SLALP,SWALP,VPR,GAM00R;
  double SLGAM00,SL2GAM00,ETAGAMR,M0R,STM0,SLM0,ETAMR,ZET1R,ETAZET,SLZET1;
  double VSBTR,SLVSBT,A1R,STA1,SLA1,SWA1,A2R,SLA2,SWA2,A3R;
  double SLA3,SWA3,COL,NTR,NFR,NFAR,NFBR,NFCR;
  double l,w,DTA;

	// Parameter information
	static ParmInfo pinfo[];

	//hyp functions
	AD hyp1(AD,double);
	AD hyp1(AD,AD);
	AD diff_hyp1(AD,double);
	AD hyp2(AD,AD, double );
	AD diff_hyp2(AD,AD, double);
	AD hyp3(AD,AD,double);
	AD diff_hyp3(AD,AD,double);
	AD hyp4(AD,double, double);
	AD diff_hyp4(AD,double, double);
	AD hyp5(AD,AD, AD);
	AD diff_hyp5(AD,AD,double);
};

#endif
