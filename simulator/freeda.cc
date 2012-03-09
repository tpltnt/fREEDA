// ----------------------------------------------------------------------
// Main routine for fREEDA
// ----------------------------------------------------------------------

extern "C"
{
  #include <stdlib.h>
  #include <stdio.h>
  #include <time.h>
  #include "./inout/parser.h"
  #include "./inout/report.h"
}

#include <iostream>
#include <cstring>
#include "./network/CircuitManager.h"
#include "help.h"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

void An_DoEndOutput(); // Defined in an_output.c
void anMain(Circuit* cir);  // Defined in anMain.cc

// Global pointer to version string
const char* freeda_version="2.0 ";

// Global pointer to the main circuit
Circuit* main_cir;

int maimed = 0;

int main(int argc, char **argv)
{
  #ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  #endif

  processCommandLine(argc, argv);

  // Open input/output files
  Init();

  // Create main circuit
  main_cir = the_CM->getCircuit("Main");

  // Parse the input netlist
  std::cout << "\n\n   *** Parsing input netlist ...\n\n";
  if (yyparse() || maimed)
  {
    report(FATAL, "Errors in netlist prevent continuation.\n\n");
  }
  try
  {
    // Print information about netlist if requested. (before expansion)
    dumpOptions();

    // Expand the main circuit (and init elements in expanded subcircuits).
    std::cout << "   *** Expanding subcircuits ...";
    fprintf(output_F, " *** Expanding subcircuits ...");
    main_cir->expand();
    main_cir->check();
    std::cout << " done.\n\n";

    // Initialize all circuits
    std::cout << "   *** Initializing and Expanding Elements ...";
    fprintf(output_F, " *** Initializing and Expanding Elements ...");
    main_cir->init();
    std::cout <<" done.\n\n";

    // Print information about netlist if requested. (output expansion)
    dumpOptions();

    // Check that reference terminals are correct
    std::cout << "   *** Checking reference terminals ...";
    fprintf(output_F, " *** Checking reference terminals ...");
    main_cir->checkReferences();
    std::cout << " done.\n\n";
    fprintf(output_F, " done.\n\n");
  }
  catch (string& error)
	{
    report(FATAL, error.c_str());
  }

  // Transfer control to the main analysis routine
  anMain(main_cir);

  // Delete the main circuit (the only one after the expansion).
  the_CM->deleteCircuit("Main");

  // Print end time
  {
    time_t curtime;
    char str[90];
    time(&curtime);
    strftime(str, 80, "%c", localtime(&curtime));
    std::cout << "\n********** fREEDA " << freeda_version << "stopping on "
              << str << " **********\n",
    fprintf(output_F, "\n\n********** fREEDA %s stopping on %s **********\n",
            freeda_version, str);
  }

  // Clean all the symbol tables, sweeps and expressions.
  majorCleanUp();

  #ifdef HAVE_MPI
  MPI_Finalize();
  #endif

  exit(0);
}

