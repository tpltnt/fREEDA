//  Nand Latch
//
//               VDD 1
//               o
//               |
//               |
//         --------------
//   -     |            |
// 2 S o---|            |---o Q 4
//   -     | Nand Latch |     - 
// 3 R o---|            |---o Q 5
//         --------------
//                |
//                |
//                o GND 6



#ifndef CmosLatchSRnand_h
#define CmosLatchSRnand_h 1

class CmosLatchSRnand : public Element
{
	public:
	       // This is the prototype of the creator routine
  		CmosLatchSRnand(const string& iname);

		// This is the prototype of the destructor routine
  		~CmosLatchSRnand() {}
         static const char* getNetlistName()
	 {
		return einfo.name;
	 }

          // local initialization
	virtual void init() throw(string&);
            private:

               // Set up the variable to store element information
		  static ItemInfo einfo;
		// Set up the variable to store the number of parameters of this element
		   static const unsigned n_par;

                  // Parameter Variables
		   double  ln, wn, lp, wp;

                 // Variables
			double C,Rs;
       
		//Parameter information. Space is allocated for the pointer to the pinfo vector      
        	static ParmInfo pinfo[];
};
#endif
