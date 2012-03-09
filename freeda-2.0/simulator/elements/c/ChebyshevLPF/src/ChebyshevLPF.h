// This is a 3 terminal highpass or lowpass Chebyshev filter model.
// It is considered as a 3 terminal device with 1 local reference node.
//              
//                 |----|        |----|
//       1 o----+--|    |-----+--|    |--------o 3
//              |  |----|     |  |----|    
//              |   g2        |    g4        
//           g1---         g3---
//             ---           ---          
//              |             |
//              |             |          
//       2 o----+-------------+                  
//         |                                   
//         ^ Local reference terminal       
//
// by Michael Bucher

#ifndef ChebyshevLPF_h
#define ChebyshevLPF_h 1

class ChebyshevLPF : public Element
{
public:
  
  ChebyshevLPF(const string& iname);

  ~ChebyshevLPF() {}

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
  // t0 t1 t2
  //
  //     ^     
  //     LRN1 
  //
  // local_ref_vec contains: {1}
  //
  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec,
			      TerminalVector& term_list);

  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);
  virtual void fillMNAM(TimeMNAM* mnam);

private:

  // Work variables
  double olda,oldb,oldg,alpha,gamma,b,g,C,L;


  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double f0, z0, ripple;
  int n, type;
  

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
