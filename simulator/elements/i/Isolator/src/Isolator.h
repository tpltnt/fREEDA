// This is an ideal analog isolator model.
// (in future variations a microwave model might be implemented)
//
//                   iso
//                ---------
//       1 o-----|  -> OK! |-----o 2
//         z->   |         |  <-z
//       3 o-----|  NO! <- |-----o 4
//                ---------
// by Don Widdows, Isac Lima, Daryl Lindsey

#ifndef Isolator_h
#define Isolator_h 1

class Isolator : public Element
{
	public:

  Isolator(const string& iname);
  ~Isolator() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // This element would add extra equations to the MNAM if
  // the reactive component to the device's impedance was implemented
  /*virtual unsigned getExtraRC(const unsigned& eqn_number,
	const MNAMType& type);
  virtual void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const;
  */

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

  //may be useful for a future microwave isolator model
  //unsigned my_row;
  //IntVector my_row;

  // Parameter variables
  double r,g;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif

