// This is Center Tap Transformer Model [3-ports]
// with 2 local reference nodes (LRN)
// 
//		I1		I3
//   n1 o->-- <=k12=>	-->----o n3 OUT1-in phase
//		|		|  		N2
//	    o	C		C  o	
//	N1	C		C	I4
//		C		C---o n4 LRN2
//		C		C
//  		C		C o 		N3
//	I2	|		| 	I5
//   n2 o----  <=k13=>	-------o n5 OUT2-out of phase
//	LRN1
//
// Olga Andreescu

#ifndef TransformerCT_h
#define TransformerCT_h 1

class TransformerCT : public Element
{
public:

TransformerCT(const string& iname);
~TransformerCT() {}

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
// Terminals ->	0	1	2	3	4
// local_ref_vector contains {1, 3}
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
DenseIntVector my_row;

// Required Parameter variables
// n1 = number of turns on primary winding
// n2 = number of turns on 1st second winding
// n3 = number of turns on 2nd second winding
// R = leakage resistance
double n1, n2, n3, R;
double k12, k13;

// Parameter information
static ParmInfo pinfo[];

};

#endif

