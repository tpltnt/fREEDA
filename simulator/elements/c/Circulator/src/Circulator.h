// This is a 3-port circulator model:
//
//
//                           / 2
//                          /
//                _________/
//              |   ->    ) |
//       1 -----|  (    <-  |
//              \ _________
//                           \ 3
//
//
//  It is essentially 3 gyrators in series with the terminals
//  labeled as such:
//
//
//                  gyr
//         Pt 1  ---------
//  t1  1 o-----|  -> pi  |-----o 2 --------
//         z->  |         |  <-z           |
//  t4  3 o-----|   0 <-  |-----o 4 -----   |  t2
//               ---------            t5 |  |
//   ------------------------------------   |
//  |  -------------------------------------
//  | |              gyr
//  | |     Pt 3  ---------
//  |  - 1 o-----|  -> pi  |-----o 2 -------
//  |       z->  |         |  <-z          |  t3
//   --- 3 o-----|   0 <-  |-----o 4 ----   |
//                ---------          t6  |  |        { t=circulator terminal }
//   ------------------------------------   |
//  |  -------------------------------------
//  | |
//  | |              gyr
//  | |     Pt 2  ---------
//  |  - 1 o-----|  -> pi  |-----o 2
//  |       z->  |         |  <-z      -> to input of first gyrator (t1, t4)
//   --- 3 o-----|   0 <-  |-----o 4
//                ---------


#ifndef Circulator_h
#define Circulator_h 1

class Circulator : public Element
{
	public:

  Circulator(const string& iname);

  ~Circulator() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // This element would add equations to the MNAM
  // if a microwave circulator was implemented (including reactances)

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
  double r,g;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif

