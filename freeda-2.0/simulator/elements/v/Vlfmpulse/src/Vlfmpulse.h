// This may look like C code, but it is really -*- C++ -*-
//
// This is a linear frequency modulation chirp source 
// (transient analysis only)
//
//       n1     + ---  -   n2
//          o----(   )-----o
//                ---     
// This element behaves as a short circuit for AC/HB analysis.
//
//
// by Frank P. Hart 5/19/2004

#ifndef Vlfmpulse_h
#define Vlfmpulse_h 1

class Vlfmpulse : public Element
{
  public:

    Vlfmpulse(const string& iname);

    ~Vlfmpulse() {}

    static const char* getNetlistName()
    {
      return einfo.name;
    }

    // This element adds equations to the MNAM
    virtual unsigned getExtraRC(const unsigned& eqn_number, 
        const MNAMType& type);
    virtual void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const;

    // fill MNAM
    virtual void fillMNAM(FreqMNAM* mnam);
    virtual void fillMNAM(TimeMNAM* mnam);
    virtual void fillSourceV(TimeMNAM* mnam);

    // State variable transient analysis
    virtual void svTran(TimeDomainSV *tdsv);
    virtual void deriv_svTran(TimeDomainSV *tdsv);

  private:

    // row assigned to this instance by the FreqMNAM
    unsigned my_row;

    // Time for the integer number of periods before current time
    double int_per_time;

    // Element information
    static ItemInfo einfo;

    // Number of parameters of this element
    static const unsigned n_par;

    // Parameter variables
    double vo, va, td, fo, deltaf, phi, tau, per;
    int chirpdir; 

    // Parameter information
    static ParmInfo pinfo[];
};

#endif

