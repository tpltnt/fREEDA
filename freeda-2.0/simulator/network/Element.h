// Element is the base class for all the possible elements.
// by Carlos E. Christoffersen

#ifndef Element_h
#define Element_h 1

#include "GraphNode.h"
#include "Terminal.h"
#include "ElementData.h"
#include <cassert>

class Circuit;
class FreqMNAM;
class TimeMNAM;
class TimeDomainSV;
class FreqDomainSV;
class WaveletDomainSV;

// ----------------------------------------------------------------------
// Flag definitions
// ----------------------------------------------------------------------
//
// Flags are used to classify elements. For example, one can request all
// the elements with the flag "linear" set to 1.
// I think the method used to implement them is quite portable and very
// efficient. The main drawback is that the number of flags is limited
// by the size of the long unsigned (at least 32 for modern computers).

typedef long unsigned ElemFlag;

// Each flag uses a bit of the ElemFlag. By now, there are only a few flags.

enum
{
  LINEAR = 1<<0,
  NONLINEAR = 1<<1,
  ONE_REF = 1<<2,
  MULTI_REF = 1<<3,
  TR_TIME_DOMAIN = 1<<4,
  TR_FREQ_DOMAIN = 1<<5,
  SOURCE = 1<<6,
  CONVOLUTION = 1<<7
};

enum MNAMType {FREQ_DOMAIN, TIME_DOMAIN};

//-----------------------------------------------------------------------

class Element : public GraphNode, public NetListItem
{
  public:

  Element(ItemInfo* einfo, ParmInfo* param_desc, const int& numparms,
  const string& iname);

  virtual ~Element();

  //---------------------------------------
  // Flags
  //---------------------------------------
  // Return true if the element satisfies the mask flags
  // i.e.: the flags set in the mask are also set in the
  // in the element flags (myflags).
  bool satisfies(const ElemFlag& mask) const;

  // ---------------------------------------
  // Terminal stuff
  // ---------------------------------------
  // Return the number of terminals specified by derived class.
  // It is required to have this variable set after init().
  inline unsigned getNumTerms() const
  {
    if (numterms)
      return numterms;
    else
      return getGraphNodeCount();
  }

  // Connect to terminals (and terminal to element).
  void connect(Terminal* terminal);

  // Get terminal by index.
  inline Terminal* getTerminal(const unsigned& i)
  {
    // We know that only Terminals can be connected to Elements.
    return (Terminal*)(getGNode(i));
  }

  // Set Circuit containing the element
  inline void setCircuit(Circuit* cir)
  {
    assert(!my_circuit);
    my_circuit = cir;
  }

  // Get the Circuit pointer
  inline Circuit* getCircuit()
  {
    assert(my_circuit);
    return my_circuit;
  }

  // Init element
  virtual void init() throw(string&);

  // Check if element has been built correctly
  virtual void check() throw(string&);

  // Get element data pointer
  inline ElementData* getElemData()
  {
    return elemdata;
  }

  // Get a vector with the indexes of the local reference nodes.
  // Get terminal pointers in term_list
  // ordered by local reference node:
  //
  // t0 t1 t2 t3 t4
  //
  //     ^        ^
  //     LRN1     LRN2
  //
  // local_ref_vec contains: {1, 4}
  //
  // This method ONLY need to be reimplemented if the local reference
  // node is not the last terminal in the list or the element has
  // multiple local reference nodes.
  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec,
  TerminalVector& term_list);


  // Propagate local external reference group to terminals belonging
  // to the same internal local reference group.
  void setGroup(unsigned ref_id);

  // Ignore setGroup()
  inline void noPropagate()
  {
    nopropagate = true;
  }

  // ---------------------------------------
  // FreqMNAM methods
  // ---------------------------------------
  // The following methods are required by the MNAM classes
  // to build the matrices.

  // Return the number of extra equations that the element requires
  // Also store internally the starting equation number for the
  // particular element instance.  Default is no extra equations, so
  // no storing is needed.
  virtual unsigned getExtraRC(const unsigned& eqn_number,
  const MNAMType& type);
  // Return starting equation number and number of rows.
  virtual void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const;

  // Fill the MNAM elements.
  virtual void fillMNAM(FreqMNAM* mnam);

  virtual void fillMNAM(TimeMNAM* mnam);
  // Elements with SOURCE flag
  virtual void fillSourceV(TimeMNAM* mnam);
  // Elements with CONVOLUTION flag
  virtual void setLastResult(DenseDoubleVector& res, const double& time);
  // Only for the Yee element
  virtual void getTimeInfo(double & step, double & stop);

  // ---------------------------------------
  // State variable related methods
  // ---------------------------------------

  // Get the number of states aported by element.
  inline unsigned getNumberOfStates() const
  {
    return my_nstates;
  }

  // Get the number of secondary state variables
  virtual unsigned getNumberOfSecStates() const;

  // ---------------------------------------
  // Analysis Methods
  // ---------------------------------------

  // State variable HB.
  int nx, ndx;
  virtual void svHB(FreqDomainSV* fdsv);
  virtual void deriv_svHB(FreqDomainSV* fdsv);

  // State variable transient.
  virtual void svTran(TimeDomainSV* tdsv);
  virtual void deriv_svTran(TimeDomainSV* tdsv);

  // State variable wavelet transient
  virtual void svWav(WaveletDomainSV* wdsv);
  virtual void deriv_svWav(WaveletDomainSV* wdsv);

  // Return delay vector for a nonlinear element
  virtual DenseDoubleVector getDelayVec();

  // ---------------------------------------
  // Methods to set and retrieve the model name
  // (to be used by the old parser code).
  // ---------------------------------------
  void setModelName(const char* name)
  {
    model_name = name;
  }

  const char* getModelName() const
  {
    return model_name.c_str();
  }

  protected:

  // Set number of state variables
  void setNumberOfStates(const unsigned& nstates);

  // Check element connections.
  virtual bool checkConnect() const;

  //---------------------------------------
  // Terminals
  //---------------------------------------
  // Return the number of terminals already connected.
  inline unsigned getTermCount() const
  {
    return getGraphNodeCount();
  }

  // Set number of terminals and allocate memory for output
  // vectors.
  void setNumTerms(const unsigned& numterms);

  //--------------------------------------------------------
  // Flags. Allow construction flags to be changed
  //--------------------------------------------------------
  inline void setFlags(const ElemFlag& flags)
  {
    myflags = flags;
  }

  private:

  // Circuit containing the element
  Circuit* my_circuit;

  // Flag to indicate that the element is being visited
  // (during graph algorithm).
  bool black;
  // Flag to avoid propagation of local reference by this element
  // To be used by elements created by the expansion of another element's
  // model.
  bool nopropagate;

  // Number of states
  unsigned my_nstates;

  //-------------------------------
  // Output variables
  //-------------------------------
  ElementData* elemdata;

  //--------------------------------------------------------
  // Flags.
  //--------------------------------------------------------
  ElemFlag myflags;

  //---------------------------------------
  // Terminals
  //---------------------------------------

  // Number of terminals.
  unsigned numterms;

  // Model name. This is a "temporary" way to support models using the
  // old existing parser code.
  string model_name;

};

// Define the element vector type.
typedef std::vector<Element*> ElementVector;

#endif
