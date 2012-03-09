// This is a v-to-v transducer implementing a Z-domain Chebyshev Type I
// lowpass or bandpass discrete-time filter.  The sampling frequency
// of the filter is fixed to the time-step, tstep, of transient
// simulation. Note that this filter will not work with variable
// time-stepping!  Also, the 'DC gain' of the filter is 1.
//
// The order of the filter is not fixed (as in AbmButterbpf10), but
// determined by the netlist specifications.  Implementation of the
// filter is in cascade form of conjugate pole pair blocks for even
// order filters.  Odd order filters put the single real pole block
// at the front of the cascade of pole pair blocks for best numerical
// properties.
//
//                 1                        3
//                  o-----+----------+-----o
//                        |          |
//           	    +  	  | Z-Domain |     +
//           	 Vin	  | Filter   |      Vout
//           	    -	  |          |     -
//                        |          |
//           	    o-----+----------+-----o
//                 2                        4
//
// by Frank P. Hart

#ifndef Zchebyi_h
#define Zchebyi_h 1

extern "C" {
#undef _C
#undef _S
}

class Zchebyi : public Element
{
  public:

    Zchebyi(const string& iname);

    ~Zchebyi();

    static const char* getNetlistName()
    {
      return einfo.name;
    }

    // Do some location initialization
    virtual void init() throw(string&) ;

    // This is a multi-reference element.
    virtual void getLocalRefIdx(UnsignedVector& local_ref_vec,
        TerminalVector& term_list);

    // State variable transient analysis.
    virtual void svTran(TimeDomainSV* tdsv);
    virtual void deriv_svTran(TimeDomainSV* tdsv);

  private:

    // Element information
    static ItemInfo einfo;

    // Number of parameters of this element
    static const unsigned n_par;

    // Parameter variables
    double fo,pbw,fls,flp,fhp,fhs,pbfdb,sbadb,ildb;
    int pord;
    bool rep;

    // Parameter information
    static ParmInfo pinfo[];

    // Passband derived parameters
    double sbw;

    // Prototype Lowpass Filter variables
    bool lowpass, Nodd, Ngiven;
    double tee, epsilon, lambda, wls, wlp, whs, whp, wp, ws, wc,
           wsh, wsl, wo2, wo;
    int N, k, kmax;
    double *plp_den1, *plp_den2, *plp_den3;
    double plp_num;
    double iltf; // insertion loss transfer function ratio

    // Z-filter variables
    double z_const, a, b, c;
    double *z_num1, *z_num2, *z_num3, *z_den1, *z_den2, *z_den3,
           *z_den4, *z_den5, *w_delp, *w_del1, *w_del2, *w_del3, *w_del4;

    // Evaluation variables
    bool do_filter_init;		// init filter one time only...
    bool first_iteration;         // update regs only at start of Newt it
    double x[2], vp[2], ip[2];	// args to eval, set in svTran...
    double *w, *x_int, *y_int; 	// for use in eval, size varies...
    double x_last, y_last;	// for use in derivative computation

    // private functions
    void filter_init();	// Initialize the Z-domain filter coefficients/delays
    double acosh(double); // Inverse sinh and cosh functions
    double asinh(double);

    // State variable transient analysis.
    virtual void eval(double*, double*,
        double*, const TimeDomainSV*) ;

    // Filter design output cout/report()
    // const string *name_str;
    char msg[80];
    int name_len;
    char *name_ch;

    // Debugging variables
    // double jacval;
};

#endif
