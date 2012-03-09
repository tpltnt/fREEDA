#ifndef ThermalBlockRC_h
#define ThermalBlockRC_h 1

class ThermalBlockRC : public Element
{
 public:

  ThermalBlockRC(const string &iname);

  ~ThermalBlockRC() {}

  static const char* getNetlistName()
    {
      return einfo.name;
    }

  //fill MNAM
  virtual void fillMNAM (FreqMNAM* mnam);
  virtual void fillMNAM (TimeMNAM* mnam);

 private:

  //Element information
  static ItemInfo einfo;
  //Number of parameters of this element
  static const unsigned n_par;

  //Parameter variable
double dn,ds,de,dw, dt, db, kmx,kmy,kmz,kild,kbulk,lx,ly,habove,hbelow,dmx, dmy, rho,cbulk, dfactor;
  static ParmInfo pinfo[];
};

#endif
