// CMOS NAND Gate
//
//
//                                           Vdd 1
//                                              o
//                                              |
//                      *-----------------------*
//                      |                       |
//                  |---+                   |---+
//                  |                       |
// Input1 2 o--|---O|    Input2 3 o--|-----O|
//             |    |                |      |
//             |    |---+            |      |---+
//             |        |            |          |
//             |        |            |          |
//             |        -------------|----------*--o 4 Output
//             |                     |          |
//             |                     |          |
//             |                     |          |
//             |                     |      |---+
//             |                     |      |
//             |                     |------|
//             |                            |
//             |                            |---+
//             |                                |
//             |                                |
//             |                            |---+
//             |                            |
//             |----------------------------|
//                                          |
//                                          |---+
//                                              |
//                                              |
//                                              o
//                                           GND 5
//
//
//        Author:
//              Shivam Priyadarshi
//

#ifndef Cmos2NandX_h
#define Cmos2NandX_h 1

class Cmos2NandX : public Element
{
	public:
	       // This is the prototype of the creator routine
  		Cmos2NandX(const string& iname);

		// This is the prototype of the destructor routine
  		~Cmos2NandX() {}
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
