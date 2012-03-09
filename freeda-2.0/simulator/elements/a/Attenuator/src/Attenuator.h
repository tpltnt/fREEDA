// This is an attenuator model
//by Shravanthi Ippagunta

// This file sets up the class for the Attenuator element

// The following is done in case Attenuator.h is included twice. If
#include <math.h>
#ifndef Attenuator_h
#define Attenuator_h 1

class Attenuator : public Element
{
  public:
    // This is the prototype of the creator routine
    Attenuator(const string& iname);

    // This is the prototype of the destructor routine
    ~Attenuator() {}

    static const char* getNetlistName()
    {
      return einfo.name;
    }

    // Do some local initialization
    virtual void init() throw(string&);

    ////
    // Function overloading is used to fil the MNAM
    // (Modified Nodal Admittance Matrix). There are two routines here
    // of identical names.  Which is called depends on the argument type.
    // fill the frequency domain form of MNAM

    //  For the frequency domain form (complex elements)
    virtual void fillMNAM(FreqMNAM* mnam);

    //  For the tiem domain form (real elements)
    virtual void fillMNAM(TimeMNAM* mnam);

  private:

    // Set up the variable to store element information
    static ItemInfo einfo;

    // Set up the variable to store the number of parameters of this element
    static const unsigned n_par;

    // List of parameter (there is only one for this element).
    // List each of the parameters of an element here.

    double r1,r2,r3,r4,zo,alpha;
    int numberOfTerminals;

    // Parameter information. Space is allocated for the pointer to the
    // pinfo vector.
    static ParmInfo pinfo[];
};

#endif


