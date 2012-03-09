/* This is the scanner used by the NPortT element
*/ 
%{
  struct
  {
    NPortT* thisNPortT;
    unsigned ncoeff;
    unsigned coeffcount;
    unsigned trcount;
    unsigned linecount;
    char serror[80];
  } scan_NPortT;
  
%}

%option noyywrap
%option never-interactive

BSPACE [ \t]+
DIGIT [0-9]
NUMBER \-?\+?{DIGIT}+
FLOAT_NUMBER [ \t]*{NUMBER}\.?{DIGIT}*([Ee]{NUMBER})?

%s PORT
%s SPECS
%s NCOEFFN
%s DATAN
%s NCOEFFD
%s DATAD

%%

<INITIAL>"#"{BSPACE}"port:group" {
  /* Prepare to read port groups */
  BEGIN(PORT);
}

<PORT>[ \t]*{NUMBER}+":"[ \t]*{DIGIT}+ {
  /* Take the second number from string. */
  int group = atoi(strstr(yytext,":")+1);
  /* Add one more port to definition. */
  scan_NPortT.thisNPortT->addPort(group);
}

<PORT>"fmax"({BSPACE})?"="({BSPACE})? {
  BEGIN(SPECS);
}

<SPECS>{FLOAT_NUMBER} {
  scan_NPortT.thisNPortT->fmax = atof(yytext);
  BEGIN(NCOEFFN);
  /* Return to the calling routine */
  return 0;
}

<NCOEFFD,NCOEFFN>{FLOAT_NUMBER} {
  sprintf(scan_NPortT.serror, "Extra polynomial coefficient at line %d", 
	  scan_NPortT.linecount); 
  return -1;
}

<DATAN,DATAD>"Degree: "{DIGIT}+ {
  sprintf(scan_NPortT.serror, "Missing polynomial coefficient at line %d", 
	  scan_NPortT.linecount); 
  return -1;
}
  

<NCOEFFN>"Degree: "{DIGIT}+ {
  // read degree of next polynomial
  scan_NPortT.ncoeff = atoi(strstr(yytext,":")+2) + 1;
  if (scan_NPortT.trcount == scan_NPortT.thisNPortT->ntransfer) {
    sprintf(scan_NPortT.serror, "Extra transfer function at line %d", 
	    scan_NPortT.linecount); 
    return -1;
  }
  // Set size of numerator polynomial
  scan_NPortT.thisNPortT->setNumSize(scan_NPortT.trcount, scan_NPortT.ncoeff);
  scan_NPortT.coeffcount = 0;
  BEGIN(DATAN);
}

<DATAN>{FLOAT_NUMBER} {
  // read numerator coefficient
  scan_NPortT.thisNPortT->setNumCoeff(scan_NPortT.trcount, 
				      scan_NPortT.coeffcount) = atof(yytext);
  scan_NPortT.coeffcount++;
  if (scan_NPortT.coeffcount == scan_NPortT.ncoeff) 
    // all coefficients are already read
    BEGIN(NCOEFFD);
}

<NCOEFFD>"Degree: "{DIGIT}+ {
  // read degree of next polynomial
  scan_NPortT.ncoeff = atoi(strstr(yytext,":")+2) + 1;
  // Set size of numerator polynomial
  scan_NPortT.thisNPortT->setDenSize(scan_NPortT.trcount, scan_NPortT.ncoeff);
  scan_NPortT.coeffcount = 0;
  BEGIN(DATAD);
}

<DATAD>{FLOAT_NUMBER} {
  // read numerator coefficient
  scan_NPortT.thisNPortT->setDenCoeff(scan_NPortT.trcount, 
				      scan_NPortT.coeffcount) = atof(yytext);
  scan_NPortT.coeffcount++;
  if (scan_NPortT.coeffcount == scan_NPortT.ncoeff) {
    // all coefficients are already read
    // Check that the transfer function is normalized (p0 == 1)
    if (abs(atof(yytext) - one) > 1e-15) {
      sprintf(scan_NPortT.serror, "p0 = %s at line %d (the transfer function must be normalized so p0 = 1)", yytext, scan_NPortT.linecount);
      return -1;
    }
    scan_NPortT.trcount++;
    BEGIN(NCOEFFN);
  }
}

"//".* |
"!".* |
{BSPACE} {
  /* Eat up white space and comments */
}

<<EOF>> {
  if (scan_NPortT.trcount != scan_NPortT.thisNPortT->ntransfer) {
    sprintf(scan_NPortT.serror, 
	    "Missing transfer function at line %d. Read %d of %d required.", 
	    scan_NPortT.linecount, scan_NPortT.trcount, 
	    scan_NPortT.thisNPortT->ntransfer); 
    return -1;
  }
  BEGIN(INITIAL);
  YY_FLUSH_BUFFER;
  yyterminate();
}

[^ #\t\n] {
  /* An unrecognized input has been found */
  sprintf(scan_NPortT.serror, "Unrecognized character %s at line %d", 
	  yytext, scan_NPortT.linecount);
  return -1;
}

\n {
  // Count lines
  scan_NPortT.linecount++;
}

%%


