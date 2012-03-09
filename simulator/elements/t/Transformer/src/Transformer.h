// This is an Ideal Transformer Model [2-ports]
// with 2 local reference nodes (LRN)
// 		I1	I3
//   n1 o->-- 	-->--o n3
//		|	|
//		C	C	
//	N1	C	C	N2
//  		C	C	
//	I2	|	|	I4
//   n2 o----	-----o n4
//	LRN1   		LRN2
//
// Olga Andreescu

#ifndef Transformer_h
#define Transformer_h 1

class Transformer : public Element
{
public:

Transformer(const string& iname);
~Transformer() {}

static const char* getNetlistName()
{
return einfo.name;
}

// Do some local initialization
virtual void init() throw(string&);

// add extra equations to the MNA Matrix
virtual unsigned getExtraRC(const unsigned& eqn_number,
		const MNAMType& type);
virtual void getExtraRC(unsigned& first_eqn, 
		unsigned& n_rows) const;

// Get a vector with the indexes of the local reference nodes.
// Terminals -> 0		1	2	3
// local_ref_vector contains {1,3}
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

// construct matrix of the stamp with extra row&column
unsigned my_row;

// Required Parameter variables
// n1 = number of turns on first winding
// n2 = number of turns on second winding
// R = leakage resistance
double n1, n2, R;

// Parameter information
static ParmInfo pinfo[];

};

#endif

