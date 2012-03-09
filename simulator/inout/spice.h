/****************************************************************
 *
 * spice.h: definitions for spice compatibility
 *
 *
 * Author:
 *    Michael Steer
 *
 *
 ****************************************************************/


#define SPICE_H
#define SPICE_SOURCE_SIN	1
#define SPICE_SOURCE_PULSE	2
#define SPICE_SOURCE_EXP	3
#define SPICE_SOURCE_PWL	4
#define SPICE_SOURCE_SFFM	5

#define OPTIONS_DEFAULT_DEFL	100.e-6
#define OPTIONS_DEFAULT_DEFW	100.e-6
#define OPTIONS_DEFAULT_DEFAD	0
#define OPTIONS_DEFAULT_DEFAS	0
#define OPTIONS_DEFAULT_TNOM	27
#define OPTIONS_DEFAULT_NUMDGT	4
#define OPTIONS_DEFAULT_CPTIME	1.e6
#define OPTIONS_DEFAULT_LIMPTS	201
#define OPTIONS_DEFAULT_ITL1	40
#define OPTIONS_DEFAULT_ITL2	20
#define OPTIONS_DEFAULT_ITL4	10
#define OPTIONS_DEFAULT_ITL5	5000
#define OPTIONS_DEFAULT_RELTOL	0.001
#define OPTIONS_DEFAULT_TRTOL	7.0
#define OPTIONS_DEFAULT_ABSTOL	1.e-12
#define OPTIONS_DEFAULT_CHGTOL	0.01e-12
#define OPTIONS_DEFAULT_VNTOL	1.e-6
#define OPTIONS_DEFAULT_PIVREL	1.e-13
#define OPTIONS_DEFAULT_GMIN	1.e-12

#define OPTIONS_DEFAULT_DIGFREQ 10.e9
#define OPTIONS_DEFAULT_DIGDRVF 2
#define OPTIONS_DEFAULT_DIGDRVZ 20.e3
#define OPTIONS_DEFAULT_DIGOVRDRV 3
#define OPTIONS_DEFAULT_DIGIOLVL 1
#define OPTIONS_DEFAULT_DIGMNTYMX 2
#define OPTIONS_DEFAULT_DIGMNTYSCALE 0.4
#define OPTIONS_DEFAULT_DIGTYMXSCALE 1.6
#define OPTIONS_DEFAULT_TNOM 27
#define OPTIONS_DEFAULT_WIDTH 80


/*
 * Define structure to hold Job Statistics
 */

typedef struct
{
  int NUNODS,
    NCNODS,
    NUMNOD,
    NUMEL,
    DIODES,
    BJTS,
    JFETS,
    MFETS,
    GASFETS,
    NUMTEM,
    ICVFLG,
    JTRFLG,
    JACFLG,
    INOISE,
    NOGO,
    NSTOP,
    NTTAR,
    NTTBR,
    NTTOV,
    IFILL,
    IOPS,
    PERSPA,
    NUMTTP,
    NUMRTP,
    NUMNIT,
    MEMUSE,
    MAXMEM,
    COPYKNT;
  double   READIN,
    SETUP,
    DCSWEEP,
    BIASPNT,
    MATSOL,
    ACAN,
    TRANAN,
    OUTPUT,
    LOAD,
    OVERHEAD,
    TOTAL_JOB_TIME;
  /*
   * Define additional variables that we will use alot.  We do not want to
   * check the capOptionsT_P all the time.
   */
  int	   echo;
  double	   gmin;
} spice_t;

extern	spice_t spiceTable;

#define SPICE_ANALYSIS_TRANSIENT	1
#define SPICE_ANALYSIS_BIAS_POINT	2


/*
  // This table is like the capOptionTable in the lib.solaris in the
  // ece603 locker.  If we cannot figure out how to use capOptionsTable
  // then we will use this table whose members will be set by the user */

/*
 * Intergation methods
 */
#define SPTR_BACKWARD_EULER	1 /* An order 1 method */
#define SPTR_TRAPEZOIDAL	2 /* An order 2 method */
#define SPTR_GEAR_TWO		3 /* Gear, 2nd order */


/*
 * Initilization methods
 */
#define SPTR_INIT_ZEROT	1	/* Initialize using zero initial conditions
				   That is, set all voltages, charges and
				   currents to zero at the first time step. */
#define	SPTR_INIT_NORM	2	/* Initialize using Use VTerm values left over
				   from last operation.
				   Typically from last iteration */
#define SPTR_INIT_STV0	3	/* Initialize using SV0
				   (i.e. value used previous iteration) */
#define	SPTR_INIT_TRAN	4	/* Initialize using Use SV1
				   (i.e. value used previous time step) */
#define	SPTR_INIT_PRDCT	5	/* Initialize using Extrapolate using SV2 and
				   ST1 */
#define	SPTR_INIT_INIT	6	/* Initialize using IC values in device table
				   value calculated and stored in model table
			           (i.e. bias point) */
#define	SPTR_INIT_OFF	7	/* Initialize using Vterm unless element
				   instance has OFF flag set in which case use
				   0 */



#define MaximumNumberOfStateVariables	4

/* GLOBAL int NumberOfStateVariables;   Number of state variables stored */

typedef struct spiceTran {
  double SPTR_StartTime;
  double SPTR_Time;		  /* Current time (all times in seconds) */
  double SPTR_DeltaT;		  /* Current time step */
  double SPTR_StopTime;		  /* Transient analysis stop time */
  double SPTR_PrintTimeStep;	  /* .PRINT and .PLOT step size */
  double SPTR_PrintSkipTime;	  /* Time at which to start output */
  double SPTR_MaxTimeStep;	  /* Maximum allowable time step */
  double SPTR_Tolerance;
  
  double SPTR_OldDeltaT[MaximumNumberOfStateVariables];
	
  /* Past time steps */
  double SPTR_ThermlVoltge;	  /* Thermal voltage at current temperature */
  int	 SPTR_IntOrd;		  /* Order of integration method */
  int	 SPTR_IntMethod;    	  /* Integration method */
  int	 SPTR_NoIterations;	  /* No. of iterations so far */
  int	 SPTR_NoForwardTimeSteps; /* No. of forward time steps */
  int	 SPTR_NoReverseTimeSteps; /* No. of reverse time steps */
  double SPTR_AnalysisTime;	  /* Analysis compute time in seconds */
  
  
  
} spiceTran_t;

/* extern spiceTran_t spiceTran; */


struct OptionTable{

  int IntegrationType;
  double RELTOL;
  double ABSTOL;
  double CHGTOL;
  double h;   /* //initial time step */


};

