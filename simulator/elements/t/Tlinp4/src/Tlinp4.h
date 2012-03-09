// This is a 4 terminal physical transmission line model.
// It is considered as a 4 terminal device with 2 local reference nodes.
//
//              
//       1 o----+------+        +-----+-------o 3
//              |      |        |     |  
//             +-+     |y12  y12|    +-+ 
//          y11| |    ---      ---   | |y11
//             | |   ( | )    ( | )  | | 
//             +-+    ---      ---   +-+ 
//              |      |        |     |  
//       2 o----+------+        +-----+-------o 4
//         |                                  |
//         ^    Local reference terminals     ^
//
// by Carlos E. Christoffersen, Mete Ozkar and Michael Steer

#ifndef Tlinp4_h
#define Tlinp4_h 1

class Tlinp4 : public Element
{
public:
  
  Tlinp4(const string& iname);

  ~Tlinp4() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // Get a vector with the indexes of the local reference nodes.
  // Get terminal pointers in term_list
  // ordered by local reference node:
  //
  // t0 t1 t2 t3 
  //
  //     ^     ^
  //     LRN1  LRN2
  //
  // local_ref_vec contains: {1, 3}
  //
  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec,
			      TerminalVector& term_list);

  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);
  virtual void fillMNAM(TimeMNAM* mnam);

private:

  // Work variables
  double c, l, alpha_nepers;


  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double k, alpha, z0mag, fscale, tand, length, fopt;
  int nsect;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
