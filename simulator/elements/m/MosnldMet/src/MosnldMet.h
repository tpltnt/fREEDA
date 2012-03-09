//
//  MET LDMOS
//
//         Author:
//               Jiankai Chang
//               Jason Thurston
//
//         There are Three state variables.  Vgs, Vds, and Vthm.
//         Here is the Large signal circuit model that represents the LDMOS that is modeled in this code
//
//         + Vrg -         /\Qgd    Gdg       +Vrd-
//  Gate ---/\/\-----|----/+-\-----/\/\---|---/\/\------------------------------------------------------------Drain
//     ^     Rg      |    \  /            |                                                                     ^
//     |             |     \/             |                                                                     |
//     |             |   + Vgd -     |----|-------|------------|----------|                                     |
//     |           + /\            |---| +        /\           |          |                                     |
//     |         Vgs/+ \           | | |         /+ \Qds       |          |                                     |
//     |           -\- /        Ids| | | Vds     \- /         ----       \--/                                   |
//     |             \/Qgs         |\ /|          \/           /\         \/                                    |
//     |             |             | | |          |   Rdiode_0/--\       ---- Rdiode_1                          |
//     |             \             |---| -        |            |          |                                     |
//     |         Ggs /               |            |            |          |                                     |
//     |             \               |            |            |          |                                     |
//     |             /               |            |            |          |                                     |
//     |             |--------|------|------------|------------|----------|       Thermal                       |
//     |                   +  |                                                      |                          |
//     |                      \                                                      |                          |
//     |                  Vrs /                                               |------|--------|                 |
//     |                      \                                             |---|    |        | +               |
//     |                      /                                             | | | Cth|        \                 |
//     |                   -  |                                      Itherm | | |   ---    Rth/ Vth_rise        |
//     |                      |                                             |\ /|   ---       \                 |
//     |                      |                                             | | |    |        /   State(2):     |
//     |                      |                                             |---|    |        | -  Vp[2]        |
//     |                      |                                               |------|--------|    is Vth_rise  |
//     |                      |                                                      |                          |
//     |                      |                                                    |---|                        |
//     |                      |                                                    | + |                        |
//     |                      |                                                    | - | V_tsnk                 |
//     |                      |                                                    |---|                        |
//     |                      |                                                      |                          |
//     |                      |                                                    -----                        |
//   State(0):                |                                                     ---  Gnd               State(1):
//   Vp[0] is                 |                                                      -                       Vp[1]
//    Vgs -------------->  Source <--------------------------------------------------------------------------is Vgd
//

//
//
//  MET LDMOS
//
//            Drain 2
//                  o
//                  |
//                  |
//              |---+
//              |   -----o 4 Power out
// Gate 1 o-----|
//              |   -----o 5 Thermal reference
//              |---+
//                  |
//                  |
//                  o
//           Source 3
//
//
//        Author:
//               Jiankai Chang
//               Jason Thurston
//
//        Last modified: 4/17/2003


#ifndef MosnldMet_h
#define MosnldMet_h 1

#include "../../../../network/ElementManager.h"
#include "../../../../network/ADInterface.h"

class MosnldMet : public ADInterface
{
	public:

  MosnldMet(const string& iname);

  ~MosnldMet() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec, TerminalVector& term_list);

	private:

  // Generic state variable evaluation routine
	virtual void eval(AD * x, AD * effort, AD * flow);

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double RG_0, RG_1, RS_0, RS_1, RD_0, RD_1, VTO_0, VTO_1, GAMMA, VST, BETA_0,
	BETA_1, LAMBDA, VGEXP,ALPHA,VK,DELTA,VBR_0,VBR_1,K1,K2,	M1,
	M2,M3,BR,RDIODE_0,RDIODE_1,ISR,NR,VTO_R,RTH,GGS,GGD,TAU,TNOM,
	TSNK,CGST,CDST,CGDT,CTH,KF,AF,FFE,N,ISS,CGS1,CGS2,CGS3,CGS4,
	CGS5,CGS6,CGD1,CGD2,CGD3,CGD4,CDS1,CDS2,CDS3,AREA,N_FING;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
