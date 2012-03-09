/********************************************************
* Implementation of the IBIS element for Freeda.
* - Ambrish Varma, Nikhil Mahajan & Rubina Ahmed.
*********************************************************/

#include "Ibis.h"

//check # of param
const unsigned Ibis::n_par = 2;

adouble Ibis::x0_old= 0;
int Ibis::r_f_old;

/****************************/
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

const char *KEYWORDS[] =
{
	"",
	"[Pullup]",
	"[Pulldown]",
	"[POWER_clamp]",
	"[GND_clamp]",
	"[Rising Waveform]",
	"[Rising Waveform]",
	"[Falling Waveform]",
	"[Falling Waveform]"
};

struct vi_pair
{
	double v;             // voltage
	double i;             // current
};

struct vt_pair
{
	double v;             // voltage
	double t;             // time
};

struct vi_out
{
	double v1, i1;
	double v2, i2;
	double v3, i3;
	double v4, i4;
} vi_output;

struct vt_out
{
	double t1, v1;
	double t2, v2;
	double t3, v3;
	double t4, v4;
} vt_output;

const int VEC_SZ = 256;   // max number of items in each VI/VT table.
const int BUF_SZ = 100;   // max length of each IBIS line.

struct vi_pair vi_table[VEC_SZ];  // the VI table.
struct vt_pair vt_table[VEC_SZ];  // the VT table.
char line_buf[BUF_SZ];

int vi_table_count = 0;
int vt_table_count = 0;

// Function prototypes
int Ibis_vt_table_reader(string Ibis_file_name, int vt_table_selection, adouble input_voltage, double V[], double T[]);
int Ibis_vi_table_reader(string Ibis_file_name, int vi_table_selection, adouble input_voltage, double V[], double I[]);
int Spline_coeff(int n, double x[], double f[], double b[], double c[],double d[]);
int Spline_value(int n, double x[], double f[], double b[], double c[], double d[],adouble  &t, int &interval, adouble &s);

/****************************/

// Element information
ItemInfo Ibis::einfo =
{
	"ibis",
	"IBIS Model using the IBIS Specification file",
	"with pullup,pulldown,ramp-up,and ramp-down and 2 V-t curves",
	DEFAULT_ADDRESS"behavioral",
	"2003_05_15"
};

// Parameter information
ParmInfo Ibis::pinfo[] =
{
	{"Ibis_file","Ibis File name given by user", TR_STRING,true},
	{"Vcc", "Vcc! Don't you know what Vcc is?!", TR_DOUBLE, false}
};

Ibis::Ibis(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
	// Set default parameter values
	paramvalue[0] = &Ibis_file;
	paramvalue[1] = &(Vcc = 3.3);

	// Set the number of terminals
	setNumTerms(4);

	// Set flags
	setFlags(NONLINEAR | ONE_REF | TR_TIME_DOMAIN);

	// Set number of states
	setNumberOfStates(3);
}

void Ibis::init() throw(string&)
{
	DenseIntVector var(3,0);
	var[0] = 0;   //current Vgin
	var[1] = 1;   //current ttime
	var[2] = 2;   //dummy
	initializeAD(var);
}


void Ibis::eval(AD * x, AD * effort, AD * flow)
{
	AD Iout1,Kd,Ku;
	AD Ipu1, Ipd1, Ipc1, Igc1;
	int r_f;
	AD Vth, Vthr=0.05*Vcc,Vthf=0.95*Vcc, Vr0_t_var,dummy;
	int interval, n=4;
	double b[11], c[11], d[11], Vtemp[4], Itemp[4], t_to_calc[4];

	AD ttime,t_var;

	ttime =  x[1];

	int table1 =1,
	table2 =2,
	table3 =3,
	table4 =4;

	for (int ind =0; ind<=3;ind++)
	{
		Vtemp[ind]=0;
		Itemp[ind]=0;
		t_to_calc[ind]=0;
	}

	// sets out the initial values at start of the simulation

	double r_ft;

	if  ((ttime==0) && (x[0]<Vthr))
	{
		r_f=5;
		r_f_old=7;
	}

	if  ((ttime==0) && (x[0]>Vthf))
	{
		r_f=7;
		r_f_old=5;
	}
	//checks whether the input just crossed any of the thresholds

	r_ft = r_f_old;
  if (((x0_old<=Vthr) && (x[0] >Vthr)) > 0.0)
  {
    r_ft = 5.0;
    t_call = x[1].val();
  }

  if (((x0_old>=Vthf) && (x[0] <Vthf)) > 0.0)
  {
    r_ft = 7.0;
    t_call = x[1].val();
  }

	r_f=(int)r_ft;

  if (x[1] > 0.0)
    t_var = x[1] - t_call;

	AD result = (double) Ibis_vt_table_reader(Ibis_file,r_f,t_var*1e-9,Vtemp, t_to_calc);

  if (result + 1.0 > 0.0)
  {
    dummy = Spline_coeff(n, t_to_calc, Vtemp, b, c, d);
    Spline_coeff(n, t_to_calc, Vtemp, b, c, d);
    Spline_value(4, t_to_calc,Vtemp, b, c, d, t_var, interval, Vr0_t_var);
  }

  if (result < 0.0)
    Vr0_t_var = Vtemp[0];

	result = (double) Ibis_vi_table_reader(Ibis_file,table1,Vr0_t_var, Vtemp,Itemp);
	if (result == 0)
	{
		Spline_coeff(n, Vtemp, Itemp,b, c, d);
		Spline_value(n, Vtemp, Itemp,b, c, d, Vr0_t_var, interval, Ipu1);
	}
	else
	{
		Ipu1 = Itemp[0];
	}

	result = (double) Ibis_vi_table_reader(Ibis_file,table2,Vr0_t_var, Vtemp,Itemp);
	if (result == 0)
	{
		Spline_coeff(n, Vtemp, Itemp,b, c, d);
		Spline_value(n, Vtemp, Itemp,b, c, d, Vr0_t_var, interval, Ipd1);
	}
	else
	{
		Ipd1 = Itemp[0];
	}

	result = (double) Ibis_vi_table_reader(Ibis_file,table3,Vr0_t_var, Vtemp,Itemp);
	if (result == 0)
	{
		Spline_coeff(n, Vtemp, Itemp,b, c, d);
		Spline_value(n, Vtemp, Itemp,b, c, d, Vr0_t_var, interval, Igc1);
	}
	else
	{
		Igc1 = Itemp[0];
	}

	// Determine Ipc1
	result = (double) Ibis_vi_table_reader(Ibis_file,table4,Vr0_t_var, Vtemp,Itemp);
	if (result == 0)
	{
		Spline_coeff(n, Vtemp, Itemp,b, c, d);
		Spline_value(n, Vtemp, Itemp,b, c, d, Vr0_t_var, interval, Ipc1);
	}
	else
	{
		Ipc1 = Itemp[0];
	}

	// calculate Ku and Kd (Ku +Kd = 1), Kd = vgin/Vcc = x[0]/Vcc
	Kd = x[0]/Vcc;
	Ku = 1-Kd;

	// now calculate Iout1 from equation (1) in paper
  if (x[1] > 0.0)
    Iout1 = (Ku * Ipu1) + (Kd*Ipd1) + Ipc1 + Igc1;

	// this will definitely not change because this directly affect ttime
	flow[2] = x[2];
	effort[2] = Vr0_t_var;

	flow[1] = x[1];
	effort[1] = x[1];

	flow[0] = Iout1;
	effort[0] = x[0];

	x0_old=x[0];
	r_f_old=r_f;
}

/////////////////////////////////////////////////////////////////////////
////////////////IBIS Reader subprogram///////////////////////////////////
/////////////////////////////////////////////////////////////////////////

/* This program peforms two tasks:
---------------------------------------------------------------------------
Task 1:
It takes as input an IBIS file name and a search keyword.
Based on the search keyword it writes Vi table vectors as output into a new
file.
--------------------------------------------------------------------------
Task2:
It takes as a single voltage value as input .
It searches through the output file created by task 1 and gives as output 2 successive vi
pairs above and below the voltage input.
------------------------------------------------------------------------------
It also does the same with VT tables.
-----------------------------------------------------------------------------
*/

// This function will return the current in Amps. For example, given
// a string like "100.26mA", this function will return 0.10026
float parse_current_string(string i_str)
{
	int start = i_str.find("mA");

	if (start < 0)
	{
		// There is no "mA" in this string.
		// Maybe it has "A" in it... try that.
		start = i_str.find("A");
		if (start > 0)
		{
			// yes! There is an "A". Remove it.
			i_str.replace(start, 1, "");          // Here, 1 is the length of "A"
		}
		return atof(i_str.c_str());
	}
	else
	{
		// we first remove the "mA" from the string and then,
		// we convert it to a float and divide it by 1000 and
		// return the result.
		i_str.replace(start, 2, "");              // Here, 2 is the length of "mA".
		return atof(i_str.c_str())/1000.0;
	}
}


// This function will return the voltage in Volts. For example, given
// a string like "100.26mV", this function will return 0.10026 volts.
float parse_voltage_string(string v_str)
{
	// check if there is a "uV" in the string.
	int start = v_str.find("uV");
	if (start < 0)
	{
		// There is no "uV" in this string.
		// Maybe it has "mV" in it... try that.
		start = v_str.find("mV");
		if (start < 0)
		{
			// There is no "mV" in this string.
			// Maybe it has a "V" in it... try that.
			start = v_str.find("V");
			if (start > 0)
			{
				v_str.replace(start, 1, "");      // 1 is the length of "V"
			}
			return atof(v_str.c_str());
		}
		else
		{
			// yes! There is an "mV". We first remove the "mV",
			// then convert it to a float and divide it by 1000
			// and return the result.
			v_str.replace(start, 2, "");          // 2 is the length of "mV"
			return atof(v_str.c_str())/1000.0;
		}
	}
	else
	{
		// we first remove the "uV" from the string and then,
		// we convert it to a float and divide it by 1000000 and
		// return the result.
		v_str.replace(start, 2, "");              // Here, 2 is the length of "uV".
		return atof(v_str.c_str())/1000000.0;
	}
}


// This function will return the time in nS. For example, given
// a string like "100.26nS", this function will return 100.26
float parse_time_string(string t_str)
{
	int start = t_str.find("nS");
	if (start < 0)
	{
		// There is no "nS" in this string.
		// Maybe it has "S" in it... try that.
		start = t_str.find("S");
		if (start > 0)
		{
			// yes! There is an "S". Remove it.
			t_str.replace(start, 1, "");          // Here, 1 is the length of "A"
		}
		return atof(t_str.c_str())*1000000000.0;
	}
	else
	{
		// we first remove the "nS" from the string and then,
		// we convert it to a float return the result.
		t_str.replace(start, 2, "");              // Here, 2 is the length of "nS".
		return atof(t_str.c_str());
	}
}

// set the vt_output structure with the vector indices passed in.
void set_vt_output(int i1, int i2, int i3, int i4)
{
	vt_output.t1 = vt_table[i1].t;
	vt_output.v1 = vt_table[i1].v;
	vt_output.t2 = vt_table[i2].t;
	vt_output.v2 = vt_table[i2].v;
	vt_output.t3 = vt_table[i3].t;
	vt_output.v3 = vt_table[i3].v;
	vt_output.t4 = vt_table[i4].t;
	vt_output.v4 = vt_table[i4].v;
}


void assign_VT_arrays(double V[4], double T[4])
{
	V[0]=vt_output.v1;
	V[1]=vt_output.v2;
	V[2]=vt_output.v3;
	V[3]=vt_output.v4;
	T[0]=vt_output.t1;
	T[1]=vt_output.t2;
	T[2]=vt_output.t3;
	T[3]=vt_output.t4;
}


/*
Function: Ibis_vt_table_reader

Inputs: Ibis_file_name: the name of the IBIS file (full path if necessary)
vt_table_selection: integer in the range [5,8] inclusive.
input_time: time in seconds. (NOT nanoseconds and NOT milliseconds !!!!)
*/
int Ibis_vt_table_reader(string Ibis_file_name, int vt_table_selection, AD input_time_sec,double V[4], double T[4])
{
	string IBIS_keyword;
	int keyword_len;
	float input_V_fixture;

	AD input_time = input_time_sec*1e9;      // to make it nanoseconds

	// reset counters!!
	vt_table_count = 0;

	// Validate Selection
	if (vt_table_selection <5 || vt_table_selection> 8)
	{
		// cout << "ERROR: Bad Selection!" <<endl;
		exit(0);
	}

	// depending on the selection, set the input V_fix
	switch (vt_table_selection)
	{
		case 5: input_V_fixture = 0.0; break;
		case 6: input_V_fixture = 3.3; break;
		case 7: input_V_fixture = 0.0; break;
		case 8: input_V_fixture = 3.3; break;
		default: input_V_fixture = -1; exit(1); break;
	}

	IBIS_keyword = string(KEYWORDS[vt_table_selection]);
	keyword_len = IBIS_keyword.length();
	// cout << " keyword Selected = " << IBIS_keyword <<endl ;

	// Read in IBIS file
	ifstream IBIS_fstream(Ibis_file_name.c_str());
	while (IBIS_fstream.getline(line_buf, BUF_SZ))
	{
		int result;
		string line_str(line_buf);
		// Print the line just for debugging purposes.
		// cout << IBIS_keyword << ":" << line_str << ":" << endl;

		// for each line, look and see if we have our keyword
		result = line_str.find(IBIS_keyword);
		if (result >= 0)
		{
			// We may have found our table. That depends on
			// what the V_fix value is.
			// R_fix
			IBIS_fstream.getline(line_buf, BUF_SZ);
			// V_fix
			IBIS_fstream.getline(line_buf, BUF_SZ);
			line_str = string(line_buf);
			istringstream line_in(line_str);
			string name, equal_sign;
			float value;
			line_in >> name >> equal_sign >> value;
			// next line only for debugging purposes...
			// cout << name << equal_sign << value << endl;
			if (value == input_V_fixture) break; // we found our VT table.
		}
	}

	// read the entire table.
	while (IBIS_fstream.getline(line_buf, BUF_SZ))
	{
		string v_str, t_str;
		string line_str(line_buf);                // got a string from a char*

		istringstream line_in(line_str);
		// Got an istringstream object. The istringstream
		// class is very convenient when you already have a string
		// and want to do formatted input from it. (like reading
		// float and int values from it just like you use cin).

		// if the very first char is a '|', ignore the line.
		if (line_buf[0] == '|') continue;

		// if the very first char is a '[', we have reached
		// another keyword. So we must stop.
		if (line_buf[0] == '[') break;

		// Print the line just for debugging purposes.
		//cout << line_str << ":" << endl;

		// Make a vector of VT values from the table.
		// To do this we have to read and understand the numbers
		// from the line_str (including the 'mV' and 'nS' notation!!)

		// first read in the voltage and current as strings
		line_in >> t_str >> v_str;
		// cout<<"time = "<< t_str <<" - voltage = "<< v_str << endl ;

		// If either of them are empty strings, we ignore this line.
		if (v_str == "" || t_str == "") continue;

		// enter this vi item into the table.
		vt_table[vt_table_count].t = parse_time_string(t_str);
		vt_table[vt_table_count].v = parse_voltage_string(v_str);
		vt_table_count++;
	}
	// By now, the vector of VT values is ready.

	// cout << "input_time1 to Ibis reader " <<input_time<< endl;

	// First check the bounds and validate the input time.
	// Since we need FOUR data points around the input time,
	// the input time has to lie between the first two times
	// and the last two times in the vt table.

	// We have to check several corner cases:

	// Case 1: input is less than the first item in the table.
	if (input_time < vt_table[0].t)
	{
		set_vt_output( 0,0,0,0 );
		assign_VT_arrays( V, T);
		return -1;
	}

	// Case 2: input is greater than the last item in the table.
	if (input_time > vt_table[vt_table_count-1].t)
	{
		set_vt_output( vt_table_count-1, vt_table_count-1,
		vt_table_count-1, vt_table_count-1 );
		assign_VT_arrays( V, T);
		return -1;
	}

	// Case 3: input is the first item in the table
	if (input_time == vt_table[0].t)
	{
		set_vt_output( 0, 1, 2, 3);
		assign_VT_arrays( V, T);
		return 0;
	}

	// Case 4: input is the last item in the table
	if (input_time == vt_table[vt_table_count-1].t)
	{
		set_vt_output ( vt_table_count-2, vt_table_count-1,
		vt_table_count-1, vt_table_count-1 );
		assign_VT_arrays( V, T);
		return 0;
	}

	// Case 5: input lies between first and second items in the table
	if (input_time > vt_table[0].t && input_time < vt_table[1].t)
	{
		set_vt_output ( 0, 1, 2, 3 );
		assign_VT_arrays( V, T);
		return 0;
	}

	// Case 6: input lies between the last and last-but-1 items in the table
	if (input_time > vt_table[vt_table_count-2].t && input_time < vt_table[vt_table_count-1].t)
	{
		set_vt_output ( vt_table_count-3, vt_table_count-2,
		vt_table_count-1, vt_table_count-1 );
		assign_VT_arrays( V, T);
		return 0;
	}

	for (int index = 1; index <= vt_table_count-3 ; index++)
	{
		// Case 7: An item is an exact match.
		if (vt_table[index].t == input_time)
		{
			set_vt_output ( index-1, index, index, index+1 );
			assign_VT_arrays( V, T);
			return 0;
		}

		// Case 8: input is between the current and the next item in the table.
		if (input_time > vt_table[index].t && input_time < vt_table[index+1].t)
		{
			set_vt_output ( index-1, index, index+1, index+2 );
			assign_VT_arrays( V, T);
			return 0;
		}

	}

	// Case 9: input is exactly equal to the the last-but-1 item.
	if (input_time == vt_table[vt_table_count-2].t)
	{
		set_vt_output ( vt_table_count-3, vt_table_count-2,
		vt_table_count-2, vt_table_count-1 );
		assign_VT_arrays( V, T);
		return 0;
	}
	return 0;
}


// set the vi_output structure with the vector indices passed in.
void set_vi_output(int i1, int i2, int i3, int i4)
{
	vi_output.i1 = vi_table[i1].i;
	vi_output.v1 = vi_table[i1].v;
	vi_output.i2 = vi_table[i2].i;
	vi_output.v2 = vi_table[i2].v;
	vi_output.i3 = vi_table[i3].i;
	vi_output.v3 = vi_table[i3].v;
	vi_output.i4 = vi_table[i4].i;
	vi_output.v4 = vi_table[i4].v;

}


void assign_VI_arrays(double V[4], double I[4])
{
	V[0]=vi_output.v1;
	V[1]=vi_output.v2;
	V[2]=vi_output.v3;
	V[3]=vi_output.v4;
	I[0]=vi_output.i1;
	I[1]=vi_output.i2;
	I[2]=vi_output.i3;
	I[3]=vi_output.i4;
}


/*
Function: Ibis_vi_table_reader

Inputs: Ibis_file_name: the name of the IBIS file (full path if necessary)
vi_table_selection: integer in the range [1,4] inclusive.
input_voltage: voltage in volts. (NOT milliVolts or anything else !!!!)

Return value:
zero : found a good set of of data points for bspline to work on.
-1 : input_voltage was out of range. Use I[0] as output current
and bypass bspline completely.
*/
int Ibis_vi_table_reader(string Ibis_file_name, int vi_table_selection, AD input_voltage,double V[4], double I[4])
{

	char line_buf[BUF_SZ];
	string IBIS_keyword;
	int keyword_len;

	// reset counters!!
	vi_table_count = 0;

	// Validate Selection
	if (vi_table_selection <1 || vi_table_selection> 4)
	{
		cout << "ERROR: Bad Selection!" <<endl;
		exit(0);
	}

	IBIS_keyword = string(KEYWORDS[vi_table_selection]);
	keyword_len = IBIS_keyword.length();
	// cout << " keyword Selected = " << IBIS_keyword <<endl ;

	// Read in IBIS file
	ifstream IBIS_fstream(Ibis_file_name.c_str());
	while (IBIS_fstream.getline(line_buf, BUF_SZ))
	{
		int result;
		string line_str(line_buf);
		// Print the line just for debugging purposes.
		// cout << IBIS_keyword << ":" << line_str << ":" << endl;

		// for each line, compare and see if we have our keyword
		result = IBIS_keyword.compare(0, keyword_len-1, line_str);
		if (result == 0) break;
	}

	// Now that we have found the keyword, read the table after it.
	// but watch out for comment lines!
	while (IBIS_fstream.getline(line_buf, BUF_SZ))
	{
		string v_str, i_str;
		string line_str(line_buf);                // got a string from a char*

		istringstream line_in(line_str);
		// Got an istringstream object. The istringstream
		// class is very convenient when you already have a string
		// and want to do formatted input from it. (like reading
		// float and int values from it just like you use cin).

		// if the very first char is a '|', ignore the line.
		if (line_buf[0] == '|') continue;

		// if the very first char is a '[', we have reached
		// another keyword. So we must stop.
		if (line_buf[0] == '[') break;

		// Make a vector of typical VI values from the table.
		// To do this we have to read and understand the numbers
		// from the line_str (including the 'mA' notation!!)

		// first read in the voltage and current as strings
		line_in >> v_str >> i_str;

		// If either of them are empty strings, we ignore this line.
		if (v_str == "" || i_str == "") continue;

		// enter this vi item into the table.
		vi_table[vi_table_count].v = atof(v_str.c_str());
		vi_table[vi_table_count].i = parse_current_string(i_str);
		vi_table_count++;
	}

	// First check the bounds and validate the input voltage.
	// Since we need FOUR data points around the input voltage,
	// the input voltage has to lie between the first two voltages
	// and the last two coltages in the vi table.

	// We have to check several corner cases:

	// Case 1: input is less than the first item in the table.
	if (input_voltage < vi_table[0].v)
	{
		set_vi_output( 0,0,0,0 );
		// assign_VI_arrays(V,I);
		I[0] = vi_table[0].i;
		return -1;
	}

	// Case 2: input is greater than the last item in the table.
	if (input_voltage > vi_table[vi_table_count-1].v)
	{
		set_vi_output( vi_table_count-4, vi_table_count-3, vi_table_count-2, vi_table_count-1 );
		// assign_VI_arrays(V,I);
		I[0] = vi_table[vi_table_count-1].i;
		return -1;
	}

	// Case 3: input is the first item in the table
	if (input_voltage == vi_table[0].v)
	{
		set_vi_output( 0, 0, 0, 1);
		assign_VI_arrays(V,I);
		return 0;
	}

	// Case 4: input is the last item in the table
	if (input_voltage == vi_table[vi_table_count-1].v)
	{
		set_vi_output ( vi_table_count-2, vi_table_count-1, vi_table_count-1, vi_table_count-1 );
		assign_VI_arrays(V,I);
		return 0;
	}

	// Case 5: input lies between first and second items in the table
	if (input_voltage > vi_table[0].v && input_voltage < vi_table[1].v)
	{
		set_vi_output ( 0, 1, 2, 3 );
		assign_VI_arrays(V,I);
		return 0;
	}

	// Case 6: input lies between the last and last-but-1 items in the table
	if (input_voltage > vi_table[vi_table_count-2].v && input_voltage < vi_table[vi_table_count-1].v)
	{
		set_vi_output ( vi_table_count-4, vi_table_count-3, vi_table_count-2, vi_table_count-1 );
		assign_VI_arrays(V,I);
		return 0;
	}

	for (int index = 1; index <= vi_table_count-3 ; index++)
	{
		// Case 7: An item is an exact match.
		if (vi_table[index].v == input_voltage)
		{
			set_vi_output ( index-1, index, index, index+1 );
			assign_VI_arrays(V,I);
			return 0;
		}

		// Case 8: input is between the current and the next item in the table.
		if (input_voltage > vi_table[index].v && input_voltage < vi_table[index+1].v)
		{
			set_vi_output ( index-1, index, index+1, index+2 );
			assign_VI_arrays(V,I);
			return 0;
		}

	}

	// Case 9: input is exactly equal to the the last-but-1 item.
	if (input_voltage == vi_table[vi_table_count-2].v)
	{
		set_vi_output ( vi_table_count-3, vi_table_count-2, vi_table_count-2, vi_table_count-1 );
		assign_VI_arrays(V,I);
		return 0;
	}
	return 0;
}


/////////////////////////////////////////////////////////////////////////
////////////////Spline subprogram////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// FUNCTION:     Functions for setting up and evaluating a cubic
//               interpolatory spline.
// AUTHORS:      Lawrence Shampine, Richard Allen, Steven Pruess  for
//               the text  Fundamentals of Numerical Computing
// DATE:         February 27, 1996
// LAST CHANGE:  April 3, 1998

int Spline_coeff(int n, double x[], double f[], double b[], double c[],
double d[])
{
	///////////////////////////////////////////////////////////////////////////////

	// Calculate coefficients defining a smooth cubic interpolatory spline.

	// Input parameters:
	//   n   = number of data points.
	//   x   = vector of values of the independent variable ordered
	//         so that  x[i] < x[i+1]  for all i.
	//   f   = vector of values of the dependent variable.
	// Output parameters:
	//   b   = vector of S'(x[i]) values.
	//   c   = vector of S"(x[i])/2 values.
	//   d   = vector of S'''(x[i]+)/6 values (i < n).
	//   return_value  =  0  normal return;
	//                 = -1  input n <= 1;
	//                 = -2  x vector is incorrectly ordered.
	///////////////////////////////////////////////////////////////////////////////

	//  Local variables:
	int i, k;
	double fp1, fpn, h, p;
	const double zero = 0.0, two = 2.0, three = 3.0;

	if (n <= 1) return -1;

	//  Calculate coefficients for the tri-diagonal system: store
	//  sub-diagonal in b, diagonal in d, difference quotient in c.

	b[0] = x[1] - x[0];
	if (b[0] <= zero) return -2;
	c[0] = (f[1] - f[0]) / b[0];
	if (n == 2)
	{
		b[0] = c[0];
		c[0] = zero;
		d[0] = zero;
		b[1] = b[0];
		c[1] = zero;
		return 0;
	}
	d[0] = two * b[0];
	for (i = 1; i < n-1; i++)
	{
		b[i] = x[i+1] - x[i];
		if (b[i] <= zero) return -2;
		c[i] = (f[i+1] - f[i]) / b[i];
		d[i] = two * (b[i] + b[i-1]);
	}
	d[n-1] = two * b[n-2];

	//  Calculate estimates for the end slopes.  Use polynomials
	//  interpolating data nearest the end.

	fp1 = c[0] - b[0] * (c[1] - c[0]) / (b[0] + b[1]);
	if (n > 3) fp1 = fp1 + b[0] * ((b[0] + b[1]) * (c[2] - c[1])
		/ (b[1] + b[2]) - c[1] + c[0]) / (x[3] - x[0]);
	fpn = c[n-2] + b[n-2] * (c[n-2] - c[n-3]) / (b[n-3] + b[n-2]);
	if (n > 3) fpn = fpn + b[n-2] * (c[n-2] - c[n-3] - (b[n-3] + b[n-2])

		* (c[n-3] - c[n-4]) / (b[n-3] + b[n-4])) / (x[n-1] - x[n-4]);

	//  Calculate the right-hand-side and store it in c.

	c[n-1] = three * (fpn - c[n-2]);
	for (i = n - 2; i > 0; i--)
		c[i] = three * (c[i] - c[i-1]);
	c[0] = three * (c[0]-fp1);

	//  Solve the tridiagonal system.

	for (k = 1; k < n; k++)
	{
		p = b[k-1] / d[k-1];
		d[k] = d[k] - p * b[k-1];
		c[k] = c[k] - p * c[k-1];
	}
	c[n-1] = c[n-1] / d[n-1];
	for (k = n - 2; k >= 0; k--)
		c[k] = (c[k] - b[k] * c[k+1]) / d[k];

	//  Calculate the coefficients defining the spline.

	for (i = 0; i < n-1; i++)
	{
		h = x[i+1] - x[i];
		d[i] = (c[i+1] - c[i]) / (three * h);
		b[i] = (f[i+1] - f[i]) / h - h * (c[i] + h * d[i]);
	}
	b[n-1] = b[n-2] + h * (two * c[n-2] + h * three * d[n-2]);
	return 0;
}


int Spline_value(int n, double x[], double f[], double b[], double c[],
double d[], AD &t, int &interval, AD &s)
{

	///////////////////////////////////////////////////////////////////////////////
	//  Evaluate the spline s at t using coefficients from Spline_coeff.
	//  Input parameters:
	//     n, x, f, b, c, d are defined as in Spline_coeff.
	//     t        = point where spline is to be evaluated.
	//  Output parameters:
	//     interval = index satisfying  x[interval] <= t < x[interval+1]
	//                unless t is outside data interval (see flag).
	//     s        = value of spline at t.
	//     return_value =  0  normal return;
	//                  = -1  n <= 1;
	//                  =  1  t < x[0];
	//                  =  2  t > x[n-1].
	///////////////////////////////////////////////////////////////////////////////

	//  Local variables:
	static int last_interval = 0;
	int flag, j;
	AD dt;

	if (n <= 1) return -1;
	flag = 0;

	//  Search for correct interval for t.

	if (t < x[0])
	{
		interval = 0;
		flag = 1;
	}
	if (t > x[n-1])
	{
		interval = n - 2;
		flag = 2;
	}
	if (flag == 0)
	{
		if (t >= x[last_interval])
			for (j = last_interval; j < n - 1; j++)
		{
			if (t < x[j+1])
			{
				interval = j;
				break;
			}
		}
		else
			for (j = last_interval - 1; j >= 0; j--)
		{
			if (t >= x[j])
			{
				interval = j;
				break;
			}
		}
		last_interval = interval;
	}

	//  Evaluate cubic polynomial on [x[interval] , x[interval+1]].

	dt = t - x[interval];
	s = f[interval] + dt * (b[interval] + dt * (c[interval] + dt  * d[interval]));
	return flag;
}

