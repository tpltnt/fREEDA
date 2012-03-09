/* This is the scanner used by the NPort element
*/
%{
  struct
  {
    NPort* thisNPort;
    ComplexDenseMatrix* ymat_P;
    double freq_factor;
    double real_part;
    unsigned row, col;
    unsigned linecount;
    char serror[80];
  } scan_NPort;

%}

%option noyywrap
%option never-interactive

BSPACE [ \t]+
DIGIT [0-9]
NUMBER \-?\+?{DIGIT}+
FLOAT_NUMBER [ \t]*{NUMBER}\.?{DIGIT}*([Ee]{NUMBER})?

%s PORT
%s SPECS
%s FREQ
%s ECOMPLEX

%%

<INITIAL>"#"{BSPACE}"port:group" {
  /* Prepare to read port groups */
  BEGIN(PORT);
}

<PORT>[ \t]*{NUMBER}+":"[ \t]*{DIGIT}+ {
  /* Take the second number from string. */
  int group = atoi(strstr(yytext,":")+1);
  /* Add one more port to definition. */
  scan_NPort.thisNPort->addPort(group);
}

<PORT>"#"({BSPACE})? {
  BEGIN(SPECS);
}

<SPECS>[GMk]?"hz"{BSPACE}"y"{BSPACE}"ri"{BSPACE}"r"{BSPACE}"50" {
  // Set frequency units
  char c = tolower(yytext[0]);
  switch (c) {
  case 'g':
    scan_NPort.freq_factor = 1e9;
    break;
  case 'm':
    scan_NPort.freq_factor = 1e6;
    break;
  case 'k':
    scan_NPort.freq_factor = 1e3;
    break;
  default:
    scan_NPort.freq_factor = one;
  }
  /* Prepare to read frequency */
  BEGIN(FREQ);
  /* Return to the calling routine */
  return 0;
}

<FREQ>{FLOAT_NUMBER}{BSPACE}{FLOAT_NUMBER} {
  // There is an extra element in file
  sprintf(scan_NPort.serror, "Extra matrix element at line %d",
	  scan_NPort.linecount);
  return -1;
}

<FREQ>{FLOAT_NUMBER} {
  /* Check if frequency is higher than previous */
  float freq = scan_NPort.freq_factor * atof(yytext);
  if (scan_NPort.thisNPort->vec_size &&
      freq < scan_NPort.thisNPort->freq_vec[scan_NPort.thisNPort->vec_size-1]) {
    // Wrong frequency order detected
    sprintf(scan_NPort.serror,
	    "Wrong frequency order detected at line %d", scan_NPort.linecount);
    return -1;
  }
  /* Add a new element in the frequency and matrix vectors. */
  if (!(scan_NPort.ymat_P = scan_NPort.thisNPort->newMatrix(freq))) {
    sprintf(scan_NPort.serror,
	    "Number of frequencies in file greater than max_freq");
    return -1;
  }
  scan_NPort.row = scan_NPort.col = 0;
  BEGIN(ECOMPLEX);
}

<ECOMPLEX>{FLOAT_NUMBER}{BSPACE}{FLOAT_NUMBER} {
  char* endptr = NULL;
  /* set real part */
  scan_NPort.real_part = strtod(yytext, &endptr);
  assert(endptr);
  /* Add element to matrix */
  assert(scan_NPort.row < scan_NPort.ymat_P->numCols());
  assert(scan_NPort.col < scan_NPort.ymat_P->numRows());
  (*scan_NPort.ymat_P)(scan_NPort.row,
		       scan_NPort.col) = double_complex(scan_NPort.real_part,
							atof(endptr));
/*   cout << "y(" << scan_NPort.row << ","  */
/*        << scan_NPort.col << ")= " <<
	  (*scan_NPort.ymat_P)(scan_NPort.row, scan_NPort.col) << endl;  */
  scan_NPort.col++;
  if (scan_NPort.col == scan_NPort.ymat_P->numRows()) {
    scan_NPort.col = 0;
    scan_NPort.row++;
  }
  if (scan_NPort.row == scan_NPort.ymat_P->numCols())
    BEGIN(FREQ);
}

<ECOMPLEX>{FLOAT_NUMBER} {
  /* Missing matrix element detected */
  sprintf(scan_NPort.serror, "Missing matrix element detected at line %d",
	  scan_NPort.linecount);
  return -1;
}

"//".* |
"!".* |
{BSPACE} {
  /* Eat up white space and comments */
}

<<EOF>> {
  BEGIN(INITIAL);
  YY_FLUSH_BUFFER;
  yyterminate();
}

[^ #\t\n] {
  /* An unrecognized input has been found */
  sprintf(scan_NPort.serror, "Unrecognized character %s at line %d",
	  yytext, scan_NPort.linecount);
  return -1;
}

\n {
  // Count lines
  scan_NPort.linecount++;
}

%%


