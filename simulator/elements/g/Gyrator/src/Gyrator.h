// This is a gyrator model, both analog and microwave.
// (For analog analysis, do not include reactance or susceptance
//  on parameter line.)
//
//                  gyr
//               ---------
//      1 o-----|  -> pi  |-----o 2
//        z1->  |         |  <-z2
//      3 o-----|   0 <-  |-----o 4
//               ---------
// by Don Widdows, Isac Lima, Daryl Lindsey

#ifndef Gyrator_h
#define Gyrator_h 1

class Gyrator : public Element
{
	public:

  Gyrator(const string& iname);
  ~Gyrator() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // This element adds equations to the MNAM
  virtual unsigned getExtraRC(const unsigned& eqn_number,
	const MNAMType& type);
  virtual void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const;

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

  //unsigned my_row;
  DenseIntVector my_row;

  double factor;

  // Parameter variables
  double r1,r2,l1,l2,x1,x2,f;

  int parameter_type;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif

