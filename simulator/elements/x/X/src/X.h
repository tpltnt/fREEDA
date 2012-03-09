// This is a subcircuit instance element
// by Carlos E. Christoffersen

#ifndef X_h
#define X_h 1

class X : public Element
{
public:

  X(const string& iname);

  ~X();

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // Attach the circuit definition
  void attachDefinition(Circuit *circuit);

  // Expand the subcircuit into the target circuit
  void expandToCircuit(Circuit* target_c);

private:

  // Overload this method because the number of connections is not constant.
  virtual bool checkConnect() const;

  // Pointer to subcircuit definition
  Circuit* circuit_def;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  // by now there are no parameters
  static const unsigned n_par;

  // Parameter information
  // static ParmInfo pinfo[];

};

#endif

