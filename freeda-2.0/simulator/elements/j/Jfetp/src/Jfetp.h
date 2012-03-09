#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

 #ifndef Jfetp_h
 #define Jfetp_h 1

 class Jfetp : public ADInterface
 {
	 public :
   static const char* getNetlistName()
   {
     return einfo.name;
   }

   Jfetp(const string& iname);
   ~Jfetp() {}
   virtual void init() throw(string&);

	 private:
   virtual void eval(AD * x, AD * effort, AD * flow);
   // some constants
   double vth;
   // parameter information
   static ItemInfo einfo;
   // number of parameters of this element
   static const unsigned n_par;
   // parameters used
   double af, area, beta, cgs0, cgd0, eg, fc, is, kf, lambda, pb,
   rd, rs, vt0, betatc, m, tnom, vt0tc, b, T;
   // parameter information
   static ParmInfo pinfo[];

 };

 #endif

