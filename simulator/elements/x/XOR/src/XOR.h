// XOR Gate 
//               VDD 1
//               o
//               | 
//               |
//         --------------
//         |            |
//2 In1o---|            |---o Q 4
//         |     XOR    |
//3 In2o---|            |
//         --------------
//                |
//                |
//                o GND 5


#ifndef XOR_h
#define XOR_h 1

class XOR : public Element
{
        public:
               // This is the prototype of the creator routine
                XOR(const string& iname);

                // This is the prototype of the destructor routine
                ~XOR() {}
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
                   double  ln, wn, lp, wp,cshunt;


                //Parameter information. Space is allocated for the pointer to the pinfo vector      
                static ParmInfo pinfo[];
};
#endif

