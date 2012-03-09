#ifndef help_h
#define help_h 1

// Print html catalog
void pCatalog();
void pCatalogElement(char *elementName);
void pCatalogFull();

// Print a help message 
void pHelpMessage();

// Print version and compilation date
void pVersion();

// Print licence information
void pLicence();

// Print verbose data about options, output requests and network if
// option opts set.
void dumpOptions();

#endif

