#ifndef Lspiral_h
#define Lspiral_h 1

class Lspiral : public Element
{
public:
  
  Lspiral(const string& iname);

  ~Lspiral() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);
  virtual void fillMNAM(TimeMNAM* mnam);

private:
  
  // Helper function for adding resistors, caps, and inductors
  Element* addLRC(int type, double* value, 
		  unsigned int t1, unsigned int t2) throw(string&);
  // Helper function for coupling inductors
  void coupleInductors(Element* l1, Element* l2) throw(string&);

  // Work variables
  unsigned int sid; // unique id for adding constituent elements

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double rdc, ldc, rs1, ls1, ms1, cw;
  double cox1, csub11, rsub11;
  double cox2, csub21, rsub21;
  double cox3, csub31, rsub31;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
