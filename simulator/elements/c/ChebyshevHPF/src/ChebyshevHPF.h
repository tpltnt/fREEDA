// This is a 3 terminal highpass Chebyshev filter model.
// It is considered as a 3 terminal device with 1 local reference node.
//
// Type II Cauer Highpass filter:
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
// For HPF g1, g3, g5, ... are inductors
//         g2, g4, g6, ... are capacitors
//
// by Michael Bucher April 21, 2008
//    Michael Steer  May 11, 2008

#ifndef ChebyshevHPF_h
#define ChebyshevHPF_h 1

class ChebyshevHPF : public Element
{
public:
  
  ChebyshevHPF(const string& iname);

  ~ChebyshevHPF() {}

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

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double f0, z0, ripple, ql, qc;
  int n, type;
  

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
