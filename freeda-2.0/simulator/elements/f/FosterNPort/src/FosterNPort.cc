#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "../../../../compat/standard.h"
#include "FosterNPort.h"
#include <cstdio>
#include <cstdlib>
#include <string>

// This is for a N-terminal device defined by y = H(s) = I(s) / V(s)
//
//  		     +----------+
//      2 o----------+          +---------o 3
//  		     +          +
//      1 o----------+          +---------o 4
//  		     +          +
//      0 o----------+          +---------o 5
//  	  	|    +          +    |
//  	  	|    +          +    |
//              |    +          +    |
//              |    +          +    |
//              |    +          +    |
//              |    +          +    |
//  		|    +----------+    |
//  		          |
//  		          |
//  		          o
//  		          Nth reference terminal


// Static members
const unsigned FosterNPort::n_par = 2;

// Element information
ItemInfo FosterNPort::einfo =
{
  "fosternport",
  "Foster N port transfer function model",
  "Ramya Mohan, Michael Steer",
  DEFAULT_ADDRESS"category:multiport",
  "2008_06_21"
};

// Parameter information

ParmInfo FosterNPort::pinfo[] =
{
  {"filename",
	"File containing the data for k, p, a, b, c, d", TR_STRING, true},	
  {"ports", "Number of Ports given in the data set", TR_INT, true}
};

FosterNPort::FosterNPort(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &filename;
  paramvalue[1] = &ports; 
  
	// Set flags
  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);
}

// init() function
//For reading the data set from the file, had to write this function because
//there was a problem reading the exponenet part and this works just fine
double FosterNPort::Value(char * line)
{
  double value = 0,fraction = 0,expon = 0;
  int index =0;
  double sign;
  if(line[index]=='-')
	{
    sign = -1;
    index ++;
	}
  else
    sign = 1;
  while(line[index]!='.' && line[index] != 0)
	{
		value = value * 10;
		value  += line[index] - '0';
		index++;
	}
  if(line[index] == 0)
    return value * sign;
  index ++;
  double pos = 10;

  while(line[index] !='e' && line[index] != 0)
    {
    fraction  += (double) (line[index] - '0')/ pos;
    index ++;
    pos = pos * 10;
    }
  value += fraction;

  index--;

  // If there is an exponential handle it.
  if(line[index] =='e')
    {
    double exp_sign;
    if(line[index+1] == '+')
      exp_sign = 1;
    else
      exp_sign = -1;
    index +=2;

    while(line[index] !=0)
      {
      expon = expon * 10;
      expon  += line[index] - '0';
      index++;
      }

  value *= pow(10,expon * exp_sign);
    }

  value = value * sign;
  return value;
}

void FosterNPort::init() throw(string&)
  {

  FILE* fp;

  // Open file with pole-residual data
  if( (fp = fopen(filename.c_str(), "r")) == NULL)
    throw(getInstanceName() + ": Cannot open input file \"" +
	filename +"\".");
  
  int i, j;   // port indeces
  int k;      // pole index
  char line[MAX_STRING_LEN]; // Temporary storage for inpuit string being read
  char keyword[MAX_STRING_LEN]; // Temporary storage for keyword
  char c;     // Temporary storage for a character
  int cindex; // character index in line
  int index;  // offset index in vector
  noCpoles = new int[ports*ports]; // Permanent Storage for Complex Pole count
  noRpoles = new int[ports*ports]; // Permanent Storage for Real Pole count
  
  

  // Set the number of terminals
  setNumTerms(ports + 1);

  // ports
  while ((c=getc(fp))=='*') while (getc(fp) !='\n' ); //skip over comments
  line[(cindex=0)++] = c;
  while((line[cindex++] = getc(fp))!= '=');
  line[cindex-1] = 0;
   if(strcmp("ports",line))
     {
     throw(getInstanceName() + ": Error in ports specification in file \""
       + filename +"\".");
     }
  cindex =0;
  while((line[cindex++] = getc(fp))!= ';'); // number of ports
  line[cindex-1]=0;
  int fileports = atoi(line);

  if(fileports != ports) 
     throw(getInstanceName() +
       ": number of ports does not agree with number of ports in file \""
       + filename +"\".");


  // Get number of complex poles and ports for each element
  {
  int numberOfPoles;
  int numberOfEntriesLeft = 2*ports*ports; // Must specify a real and complex
                                         // pole count for each element.
  for(i=1; i<ports; i++) // set pole counts to -1 so we can check we have them
    for(j=0; j<ports; j++)
      noCpoles[i,j] = noRpoles[i,j] = -1;

  while(numberOfEntriesLeft)
    {
    while (getc(fp) !='\n' ); //get past the linefeed
    while ((c=getc(fp))=='*') while (getc(fp) !='\n' ); //skip over comments

    keyword[(cindex=0)++] = c; // get keyword
    while((keyword[cindex++] = getc(fp)) != '(');
    keyword[cindex-1] = 0;
  
    cindex = 0; // get i index [i,j] th element
    while((line[cindex++] = getc(fp))!= ',');
    line[cindex-1] = 0;
    i = atoi(line)-1;

    cindex = 0; // get j index
    while((line[cindex++] = getc(fp))!= ')');
    line[cindex-1] = 0;
    j = atoi(line)-1;

    while(getc(fp)!= '=');

    cindex = 0; // get number of poles for element
    while((line[cindex++] = getc(fp))!= ';');
    line[cindex-1] = 0;
    numberOfPoles = atoi(line);

    index = i*ports + j;
    if(!strcmp("cpoles",keyword))
       noCpoles[index] = numberOfPoles;
    else if(!strcmp("rpoles",keyword))
       noRpoles[index] = numberOfPoles;
    else
      throw(getInstanceName() +
      ": Error in specification number of poles in file \"" + filename +"\".");

    numberOfEntriesLeft--;
    }

  }


  // Allocate Space for poles and residues
  cPoleRe = new double*[ports*ports];
  cPoleIm = new double*[ports*ports];
  cResidueRe = new double*[ports*ports];
  cResidueIm = new double*[ports*ports];
  rPole = new double*[ports*ports];
  rResidue = new double*[ports*ports];
  for(i=0; i<ports; i++)
    for(j=0; j<ports; j++)
      {
      index = i*ports + j;
      // Complex poles and residues
      if(noCpoles[index]>0)
        {
        cPoleRe[index] = new double[noCpoles[index]];
        cPoleIm[index] = new double[noCpoles[index]];
        cResidueRe[index] = new double[noCpoles[index]];
        cResidueIm[index] = new double[noCpoles[index]];
        }
      else
        {
        cPoleRe[index] = 0;
        cPoleIm[index] = 0;
        cResidueRe[index] = 0;
        cResidueIm[index] = 0;
        }
      for(k=0; k<noCpoles[index]; k++)
        {
        *(cPoleRe[index] +k) = 0.;
        *(cPoleIm[index] +k) = 0.;
        *(cResidueRe[index] +k) = 0.;
        *(cResidueIm[index] +k) = 0.;
        }
      // Real poles and residues
      if(noRpoles[index]>0)
        {
        rPole[index] = new double[noRpoles[index]];
        rResidue[index] = new double[noRpoles[index]];
        }
      else
        {
        rPole[index] = 0;
        rResidue[index] = 0;
        }
      for(k=0; k<noRpoles[index]; k++)
        {
        *(rPole[index] +k) = 0.;
        *(rResidue[index] +k) = 0.;
        }
      }


  // Get complex and real, poles and residues for each element
  {
  double realValue, imaginaryValue;

  // Calculate number of entries expected
  int numberOfEntriesLeft = 0;
  for(i=0; i<ports; i++)
    for(j=0; j<ports; j++)
      {
      index = i*ports + j;
      numberOfEntriesLeft +=  2*(noCpoles[index]); // One each for pole, residue
      numberOfEntriesLeft +=  2*(noRpoles[index]); // One each for pole, residue
      }

  while(numberOfEntriesLeft)
    {
    while (getc(fp) !='\n' ); //get past the linefeed
    while ((c=getc(fp))=='*') while (getc(fp) !='\n' ); //skip comments

    keyword[(cindex=0)++] = c; // get keyword
    while((keyword[cindex++] = getc(fp)) != '(');
    keyword[cindex-1] = 0;
  
    cindex = 0; // get i index [i,j] th element
    while((line[cindex++] = getc(fp))!= ',');
    line[cindex-1] = 0;
    i = atoi(line) -1;

    cindex = 0; // get j index
    while((line[cindex++] = getc(fp))!= ',');
    line[cindex-1] = 0;
    j = atoi(line) -1;

    index = i*ports + j;

    cindex = 0; // get j index
    while((line[cindex++] = getc(fp))!= ')');
    line[cindex-1] = 0;
    k = atoi(line) -1;

    while(getc(fp)!= '=');

    // Read in complex pole or residue value
    if(!strcmp("compp",keyword) || !strcmp("compr",keyword))
      {
      // get real part
      cindex = 0;
      while((line[cindex++] = getc(fp))!= 'i');
      line[cindex-2] = 0;
      realValue = Value(line);

      // The next character should be '*'
      if(getc(fp)!= '*')
        throw(getInstanceName() +
      ": Error in imaginary number specification in file \"" + filename +"\".");

      // get imaginary part
      cindex = 0;
      while((line[cindex++] = getc(fp))!= ';');
      line[cindex-1] = 0;

      imaginaryValue = Value(line);

      if(!strcmp("compp",keyword))
        {
        *(cPoleRe[index] +k) = realValue;
        *(cPoleIm[index] +k) = imaginaryValue;
        }
      else
        {
        *(cResidueRe[index] +k) = realValue;
        *(cResidueIm[index] +k) = imaginaryValue;
        }
      }

    // Read in real pole or residue value
    // Line can end in ';' or linefeed
    else if(!strcmp("realp",keyword) || !strcmp("realr",keyword))
      {
      cindex = 0; 
      while((line[cindex++] = getc(fp))!= ';');
        line[cindex-1] = 0;
      realValue = Value(line);

      if(!strcmp("realp",keyword))
        *(rPole[index] +k) = realValue;
      else
        *(rResidue[index] +k) = realValue;
      }

    // Check for error
    else
      throw(getInstanceName() +
      ": Error in specification of poles in file \"" + filename +"\".");

    numberOfEntriesLeft--;
    }

  }


  }

unsigned
   FosterNPort::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
  {
  int i, j, index;
  if (type == TIME_DOMAIN)
    {
    // Add extra rows and columns
    my_start_row = eqn_number;

    extra_rcs = 0;

    for(i=0; i<ports; i++)
      for(j=0; j<ports; j++)
        {
        index = i*ports + j;
        extra_rcs +=  3*noCpoles[index] + noRpoles[index];
        }

    }
  else
    {
    // No need for extra rows and columns in the frequency domain
    my_start_row = 0;
    extra_rcs = 0;
    //i_start_row = 0;
    }

  return extra_rcs;
  }

void FosterNPort::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
  {
  first_eqn = my_start_row;
  n_rows = extra_rcs;
 }

void FosterNPort::fillMNAM(FreqMNAM* mnam)
  {
  double_complex jw = mnam->getFreq() * double_complex(0., 2.) * pi;
  double rRe;  // Real part of residue
  double rIm;  // Imaginary part of residue
  double pRe;  // Real part of pole
  double pIm;  // Imaginary part of pole
  double A, B, C, D;
  int index;  // offset index in vector

  for(int i=0; i<ports; i++)
    for(int j=0; j<ports; j++)
      {
      index = i*ports + j;

      double_complex g = zero;

      // Handle complex poles
      for(int k=0; k< noCpoles[index]; k++)
        {
        rRe = *(cResidueRe[index] +k);
        rIm = *(cResidueIm[index] +k);
        pRe = *(cPoleRe[index] +k);
        pIm = *(cPoleIm[index] +k);

        A = 2 * rRe;
        B = -2 * (rRe*pRe + rIm*pIm);
        C = -2 * pRe;
        D = pRe*pRe + pIm*pIm;
	g += (A*jw + B)/(jw*jw + C*jw +D);
        }

      // Handle real poles
      for(int k=0; k< noRpoles[index]; k++)
        {
        rRe = *(rResidue[index] +k);
        pRe = *(rPole[index] +k);
	g += rRe / (jw - pRe);
        }

      // Incorporate in MNAM
      if(g != double_complex(0))
        {
        unsigned refIndex = getTerminal(ports)->getRC();
        unsigned iIndex = getTerminal(i)->getRC();
        unsigned jIndex = getTerminal(j)->getRC();
        mnam->setQuad(iIndex, refIndex, jIndex, refIndex ,g);
        }
      }
  }


void FosterNPort::fillMNAM(TimeMNAM* mnam)
  {
  int i, j;     // port indices
  int k;        // pole indices
  int index;    // offset index in vector
  unsigned rowcol;   // row/column index into MNA
  unsigned refIndex; // row/column index into MNA
  unsigned jIndex;   // row/column index into MNA
  unsigned iIndex;   // row/column index into MNA
  double rRe;   // Real part of residue
  double rIm;   // Imaginary part of residue
  double pRe;   // Real part of pole
  double pIm;   // Imaginary part of pole
  double A, B, C, D;

  cout << "Filling Time Domain MNAM for " << getInstanceName()  << endl;

  rowcol = my_start_row; // get index of the first extra equation
  for(int i=0; i<ports; i++)
    for(int j=0; j<ports; j++)
      {
      index = i*ports + j;

      // Add MNAM entries from real poles
      for(int k=0; k < noRpoles[index]; k++)
        {
        rRe = *(rResidue[index] +k); // There is only a real part
        pRe = *(rPole[index] +k);    // There is only a real part
        refIndex = getTerminal(ports)->getRC();// get index of reference
                                               // terminal row 
        jIndex = getTerminal(j)->getRC();// get index of j th terminal row 
        iIndex = getTerminal(i)->getRC();// get index of i th terminal row 

        // Set entries in MNAM, MNAM Derivative
          {
          // Set MNAM row
          mnam->setMElement(rowcol, jIndex, rRe);
          mnam->setMElement(rowcol, refIndex, -rRe);
  
          // Set  MNAM column
	  mnam->setMElement(iIndex, rowcol, one);
	  mnam->setMElement(refIndex , rowcol, -one);
  
	  mnam->setMElement(rowcol, rowcol, pRe);

        // Set entries in derivative of MNAM, Mp
          mnam->setMpElement(rowcol, rowcol, -one);
  
          rowcol++; // move onto next extra equation in MNAM
          }
        }


      // Add MNAM entries from complex poles
      // Need to enter two rows/columns in MNAM
      for(int k=0; k < noCpoles[index]; k++)
        {
        rRe = *(cResidueRe[index] +k); // Get the real part of the residue
        rIm = *(cResidueIm[index] +k); // Get the imaginary part of the residue
        pRe = *(cPoleRe[index] +k);    // Get the real part of the pole
        pIm = *(cPoleIm[index] +k);    // Get the imaginary part of the pole

        refIndex = getTerminal(ports)->getRC();// get index of reference
                                               // terminal row 
        jIndex = getTerminal(j)->getRC();// get index of j th terminal row 
        iIndex = getTerminal(i)->getRC();// get index of i th terminal row 

        A = 2 * rRe;
        B = -2 * (rRe*pRe + rIm*pIm);
        C = -2 * pRe;
        D = pRe*pRe + pIm*pIm;

        // Set entries in MNAM, MNAM Derivative
          {
          // Set MNAM row 1
          mnam->setMElement(rowcol, jIndex, -B);
          mnam->setMElement(rowcol, refIndex, B);
          mnam->setMElement(rowcol, rowcol, D);
          mnam->setMElement(rowcol, rowcol+1, C);
          mnam->setMElement(rowcol, rowcol+2, one);
          // Set MNAM row 2
          mnam->setMElement(rowcol+1, rowcol+1, one);
          // Set MNAM row 3
          mnam->setMElement(rowcol+2, rowcol+2, one);
          // Set  MNAM column 1
	  mnam->setMElement(iIndex, rowcol, one);
	  mnam->setMElement(refIndex , rowcol, -one);
  
          // Set Mp row 1
          mnam->setMpElement(rowcol, jIndex, -A);
          mnam->setMpElement(rowcol, refIndex, A);
          // Set Mp row 2
          mnam->setMpElement(rowcol+1, rowcol, -one);
          // Set Mp row 2
          mnam->setMpElement(rowcol+2, rowcol+1, -one);
  
          rowcol +=3; // move onto next extra equatioa setn in MNAM
          }
        }
     }
  rowcol--;
  cout << "Complete Time Domain MNAM for " << getInstanceName() << ", number of equations = " << rowcol  << endl;
  }
