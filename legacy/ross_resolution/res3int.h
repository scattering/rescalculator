/*>>>>  global input data declarations and dictionary setup  <<<<<<<<<*/

//*>>>>>>>>>>>>>> forward function declarations <<<<<<<<<<<<<*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* convolution data structs */
#include "conv3.h"

/* external entry point declarations */
#include "res3.h"

/* multixtal functions */
#include "multixtal.h"

/* for storing and recalling res in and out by var name string */
#include "dictp.h"

#include <sys/types.h>
#include <ctype.h>
#include <sys/uio.h>


#include "nonlineq.h"   //* nonlinear equs stuff for peak and width fitting */
#include "lreduce.h"    //* matrix alloc & lin equations reduction declarations */
#include "vector.h"     //* vector function declarations */


static int set_recip() ;
static int sampangles(double scan[], double as[], double bs[], double cs[],
		      double au[], double ap[], double ki, double kf,
		      double ang[]) ;

static int set_cs() ;
static int set_fix() ;

static int set_spin() ;
static int set_xaxis() ;

static int width_calc() ;
static int lowercase(char *) ;

static int recip() ;
static void transform(double TM[3][3], double *in, double *out) ;

static void lorentz( double k, double *lorfac, double *rlorfac ) ;

static int resm( double eI, double eF, double q ) ;


static int geocs( double x, double *v ) ;
static int hhgeocs( double x, double *v ) ;
static int geocshkle( double scan[], double *v ) ;
static int hhhh( double x, double *v ) ;



static double recvec( double hkl[], double uni[] ) ;

static double vrv( double v1[], double v2[] ) ;
static double vrvq( double v1[], double v2[] ) ;

static void   rmvec( double v1[], double v2[] ) ;

static int calc_resm( double s[], int flag ) ;
static int calc_resm_aux( double scan[] ) ;
static int plotstuff() ;

static int ellax2( double **m, double p, double hw[2],
		   int ix, double ***pve ) ; /* pve 3,2,2 */
static int proji(  double **m, double a[3], double ***is,double ***pj) ;
/* is 3,2,2  pj 3,2,2 */

static int jacobi( double **aa, int n, double *d, double **v, int *nrot ) ;

static double etok( double ) ; double etow( double ) ;
static int kif( int ityp, double efix, double e, double *ki, double *kf ) ;


#define	 PI		3.14159265
#define  TWOPI		(2*PI)    /*6.2831853*/
#define  ISZERO         1.e-20



/*************************** res calculations ********************/

#define  RAN		0777777

#define	 PIOVER2	(PI/2.)   /*1.5707963*/

#define  RT2PI          2.506628275
#define  RT2PINV        0.39894228
#define  DEGTORAD	0.017453293
#define	 RADTODEG	57.295778
#define  NEUTSTIFF 	2.072141789	/* hbar**2/2massneut */
#define  DNEUT		2.072141789
#define  EXPARGLIM	709.		/* biggest arg of exp < infinity */
#define  GUESS1		0		/* scan point guesses for cs min */
#define	 GUESS2		2




/*

  All relevant i/o variables are global 
  since they are used in multiple functions

  A dictionary is used for input parsing 
  to store correspondence between cmnd strings
  and data pointers as well as data type

  So link with dictp.cc

  June 2000 RWE  change all 2D and 3D matrices to pointers 
  to unify calling protocol

*/


/* static EQUATION equ[2] */
/* holds linear constraint equations as C0 + SUM Ci Xi = 0 */


/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/



/*>>>>>>>>>>>>>>>> STUFF for EIGEN calculations <<<<<<<*/

//static double tol = 1.e-5 ;
static double fmprec = 1.e-3 ;
static double eps = 1.e-8 ;
//static int    maxit = 1000 ;

/*
#define	TINY	1.0e-20 ;
*/

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);


/*
 * Macros for testing floating-point values for certain special cases:
 *
 *	IS_NAN	Test for not-a-number by comparing a value against itself
 *	IF_INF	Test for infinity by comparing against the largest floating
 *		point value.
 */

#define IS_NAN(v) ((v) != (v))

#ifdef DBL_MAX
#   define IS_INF(v) (((v) > DBL_MAX) || ((v) < -DBL_MAX))
#else
#   define IS_INF(v) 0
#endif


ENTRY *globalinput ;
ENTRY *globaloutput ;
DICT *globalinputdict ;
DICT *globaloutputdict ;

static char *errmsg = NULL ;
static int errmsgAlloc = 0 ;
static char *res3memerr = "res3 memory error" ;

//* SAMPLE */
static double latt[3] = {TWOPI, TWOPI, TWOPI} ;
static double angl[3] = {90., 90., 90.} ;
static double ahkl[3] = {1, 0, 0} ;
static double bhkl[3] = {0, 0, 1} ;


static double ahklx[3], bhklx[3], phklx[3], zhklx[3] ;
static double sampleArea = 1. ;  /* in cm^2 */
static double sampleThick = 1. ; /* in cm */
static double Fsq = 1. ;         /* in barns */
static double Fsqeff = 1. ;
static double absorbCoef = 0. ; /* in cm-1 */

//* SCAN and CROSS-SECTION */
static double hkle[4] = {1, 0, 0, 0} ;
static double Ki, Kf ;
static double hklei[4] = {1, 0, 0, 0} ;
static double dopt[4] = {1, 0, 0, 0} ;
static double step[4] = {0.05, 0, 0, 0} ;
static double point[4]= {1, 0, 0, 0} ;
static double dire1[4]= {0.1, 0, 0, 10} ;
static double dire2[4]= {0, 0, 1, 0} ;
static double dire3[4]= {0, 0, 0, 1} ;
static char   cstyp[20]   =  "Bragg_point" ;
static double muoff = 0. ;

//* SPECTROMETER */
static double hcol[4] = {40, 40, 40, 40} ;
static double vcol[4] = {2,  2,  2,  2} ;
static double mosa[3] = {40, 40, 40} ;
static double vmos[3] = {-1, -1, -1} ;
static double dsps[2] = {3.35416, 3.35416} ;
static double efix    = 14.7 ;
static double defix   = 0. ;
static double efixd ;
static double lamb    = 2.359 ;
static double rxtotalflux = 0.04 ; //* NBSR total thermal flux in n/A**2/sec */
static double epeakflux = 30 ;         //* energy at peak in flux (meV)  */
static double moneff = 1.e-4 ;         //* at 1 A-1 */
static double monArea = 1. ;           //* cm**2 */
static double mref = 1. ;
static double aref = 1. ;
//static char   efixtyp[8] =  "mono" ;
static char   aspin[32]   =  "lrl" ;
static int scattsideright = 1 ;
static char *monoprog = NULL ;
static char *analprog = NULL ;
static char *datasrc = NULL ;
static double angleA = 0. ;
static double chiA = 0. ;
static double chiB = 0. ;
static double omegaoff[2] = {0., 0.} ;   
static double tilts[2]    = {0., 0.} ;
//* xtal angle offsets in degrees */

static double xposition[5]   = {0., 0., 0., 0., 0.} ;   
//* xtal position offsets, x along incident */
static double yposition[5]   = {0., 0., 0., 0., 0.} ;
static double zposition[5]   = {0., 0., 0., 0., 0.} ;

static double thick[5]      = {0., 0., 0., 0., 0.} ;   
static double width[5]       = {0., 0., 0., 0., 0.} ;
static double height[5]      = {0., 0., 0., 0., 0.} ;
//* xtal dimensions, width is along omega */

static double distance[4]    = {1., 1., 1., 1.} ;       
//* reac-mono,mono-samp,samp-anal,anal-det dist */


static int monomode = 0 ;
static int analmode = 0 ;
static int geomode = 0 ;

/*
  April 2003
  Change global angle offsets and geom coll default (mode=0)
  direct input in angles
  instead of calulated from distances and dimensions
  these are in degrees
*/
static double gammaoff[4] = { 0., 0., 0., 0. } ;
static double deltaoff[4] = { 0., 0., 0., 0. } ;
static double geoalpha[4] = { 180., 180., 180., 180. } ;
static double geobeta[4]  = { 180., 180., 180., 180. } ;

//* MISC */
static char   xaxistyp[8]=  "ahkl" ;
static double ref1[3] = {0., 0., 0.} ;
static double ref2[3] = {0., 0., 0.} ;
static double ref3[3] = {0., 0., 0.} ;
static int    reduced = 0 ;
static int    ifixmon = 1 ;
static int    offseton = 0 ;
static int    calcATpeak = 0 ;
static int    npts = 2 ;

//* Integration */
static int Imax[4] = { 8, 8, 8, 8 } ;
static int Imin[4] = { 0, 0, 0, 0 } ;
static double Irel[4] = { 1.e-4, 1.e-4, 1.e-4, 1.e-4 } ;
static double Iabs[4] = { 1.e26, 1.e26, 1.e26, 1.e26 } ;
static double resprobmin = 1.e-4 ;

//*>>>>>>>>>>>>>>>>  output variables <<<<<<<<<<<<<*/


static double vres, fwhmz ;
static double hkleI[4], hkleF[4] ;
static double angls[4], cota, tth0 ;
static double peak[4], peakmag ;
static double fwhm[4], fwhmmag ;
static double fwhmb[3], fwhmeinc, fwhmp[3] ;
static double fwhmang[4], fwhmcot ;
static double tot, cor, relint, maxcountrate ;
static double extinctLength, vol, mon, hklecountrate, rmv ;

static double ***pr = (double***)0 ;
static double ***sc = (double***)0 ;

static double ***prp = (double***)0 ;
static double ***scp = (double***)0 ;

static double axproj[3] ;
static double axprojp[3] ;

static double **sedir = (double**)0 ;
static double **pedir = (double**)0 ;

static double **sedirp = (double**)0 ;
static double **pedirp = (double**)0 ;

static double **eigenvecs = (double**)0 ;

static double eigenvals[3] ;
static double xyeoff[3], qzoff, hkloff[3], hkleff[4] ;
/* static double xyecnoff[3], qzcnoff */
static double kixy[2], kfxy[2], kimag, kfmag ;
static double xcnoff, ycnoff, zcnoff, ecnoff ;
static double kioff, kfoff ;
static double kie, kfe, kinom, kfnom ;

static double **rm = (double**)0 ;
static double **rmcn = (double**)0 ;

static double avec[3], bvec[3], cvec[3] ;
static double zuni[3] ;
static double ast[3], bst[3], cst[3], dpstar[6] ;
static double rlat[3], rang[3] ;
static double auni[3], amag ;
static double apun[3], bmag ;
static double apub[3], abdot, apbdot ;
static double aunix[3], apunx[3], zunix[3], bunix[3] ;
static double q0[2], ref1xy[2], ref2xy[2], ref3xy[2] ;
static double ref1mag, ref2mag, ref3mag ;
static double cspoint[4], csdire1[4], csdire2[4], cs1diru[2] ;
static double hklmag, ihklmag, stepmag, pointmag, diremag[3] ;
static double hklqpun[3], qpmag ;

static int nfree ;

static double L[3][3] ; /* dotpros of recip latt vecs */
static double M[3][3] ; /* dotpros of real space latt vecs w auni apun zuni */
static double U[3][3] ; /* from chiA chiB deltaA */
static double MX[3][3] ; /* M * U */
static double MI[3][3] ; /* MXtranspose * L */

/*>>>>>>>>>>>>>>  associated output variables <<<<<<<<<<<<<<*/

static double lor, r0, cellvol ;
static double bmagp, qmag, qsqr ;
static double xmag, ymag, xuni[3], yuni[3] ;
static double buni[3] ;
static double quni[3], qpun[3], stepuni[3] ;

static double **rmp = (double**)0 ;
static double **rmw = (double**)0 ;
/* static double **rm44 = (double**)0 */
static double **rmq = (double**)0 ;
static double **rmplot = (double**)0 ;

static double xyeoffw[3] ;
static double xyeoffp[3] ;
static double detrm, rminv[3][3], brm[3] ;
static double gc[3][3], gl[3], dfh[3][3] ;
static double c2bar[3][3], bbar[3] ;

static double countrate ;

static double phkl[3], phkln[3], hkla[3], hklp[3], hklz[3] ;
//static double kisq, kfsq ;

static int    iwctyp ;

/*>>>>>>>>>>>>>>>  associated input variable <<<<<<<<<<*/

static int spin[3] = {1,1,1} ; /* lrl */
static int ifix    = 2 ; /* mono */
static int iax     = 1 ; /* ahkl */
static int icstyp  = 0 ; /* point */

/* sample fields */
static double Temp ;
static double Hfld ;
static double PolX, PolY, PolZ ;

static double *dbls = NULL ;
static int ndblalloc = 0 ;
static int incralloc = 16 ;


static double rmpM[3][3], rmwM[3][3], rmqM[3][3] ;
static double rmplotM[3][3] ;
static double rmM[3][3], rmcnM[3][3] ;
static double prM[3][2][2], scM[3][2][2] ;
static double prpM[3][2][2], scpM[3][2][2] ;
static double sedirM[2][2], pedirM[2][2] ;
static double sedirpM[2][2], pedirpM[2][2] ;
static double eigenvecsM[3][3] ;

  static DICT inputdict  = { 0, (ENTRY**)0 } ;
  static DICT outputdict = { 0, (ENTRY**)0 } ;

  /* create global array storage */

  /* input dictionary entries    name address typeflag */

  static ENTRY  inputentries[] = {
    {"LATT1",{latt},2},{"LATT2",{latt+1},2},{"LATT3",{latt+2},2},
    {"ANGL1",{angl},2},{"ANGL2",{angl+1},2},{"ANGL3",{angl+2},2},
    {"AHKL1",{ahkl},2},{"AHKL2",{ahkl+1},2},{"AHKL3",{ahkl+2},2},
    {"BHKL1",{bhkl},2},{"BHKL2",{bhkl+1},2},{"BHKL3",{bhkl+2},2},
    {"CHIA",{&chiA},2},{"CHIB",{&chiB},2},{"ANGLA",{&angleA},2},
    {"HCOL1",{hcol},2},{"HCOL2",{hcol+1},2},
    {"HCOL3",{hcol+2},2},{"HCOL4",{hcol+3},2},
    {"GHCO1",{geoalpha},2},{"GHCO2",{geoalpha+1},2},
    {"GHCO3",{geoalpha+2},2},{"GHCO4",{geoalpha+3},2},
    {"VCOL1",{vcol},2},{"VCOL2",{vcol+1},2},
    {"VCOL3",{vcol+2},2},{"VCOL4",{vcol+3},2},
    {"GVCO1",{geobeta},2},{"GVCO2",{geobeta+1},2},
    {"GVCO3",{geobeta+2},2},{"GVCO4",{geobeta+3},2},
    {"MOSAM",{mosa},2},{"MOSAS",{mosa+1},2},{"MOSAA",{mosa+2},2},
    {"VMOSM",{vmos},2},{"VMOSS",{vmos+1},2},{"VMOSA",{vmos+2},2},
    {"DSPM", {dsps},2},{"DSPA", {dsps+1},2},
    {"EFIX", {&efix},2},{"DEFIX", {&defix},2},{"LAMB", {&lamb},2},
    {"FIXM", {&ifixmon},1},
    {"SCATR",{&scattsideright},1},
    {"LR",   {aspin},0},
    {"OMOFFM",{omegaoff},2},{"OMOFFA",{omegaoff+1},2},
    {"TILTM",{tilts},2},{"TILTA",{tilts+1},2},
    {"XPOSR",{xposition},2},
    {"XPOSM",{xposition+1},2},{"XPOSS",{xposition+2},2},
    {"XPOSA",{xposition+3},2},{"XPOSD",{xposition+4},2},
    {"YPOSR",{yposition},2},
    {"YPOSM",{yposition+1},2},{"YPOSS",{yposition+2},2},
    {"YPOSA",{yposition+3},2},{"YPOSD",{yposition+4},2},
    {"ZPOSR",{zposition},2},
    {"ZPOSM",{zposition+1},2},{"ZPOSS",{zposition+2},2},
    {"ZPOSA",{zposition+3},2},{"ZPOSD",{zposition+4},2},
    {"THKR",{thick},2},{"THKM",{thick+1},2},{"THKS",{thick+2},2},
    {"THKA",{thick+3},2},{"THKD",{thick+4},2},
    {"HGTR",{height},2},{"HGTM",{height+1},2},{"HGTS",{height+2},2},
    {"HGTA",{height+3},2},{"HGTD",{height+4},2},
    {"WIDR", {width},2},{"WIDM", {width+1},2},{"WIDS", {width+2},2},
    {"WIDA", {width+3},2},{"WIDD", {width+4},2},
    {"DRM",{distance},2},{"DMS",{distance+1},2},{"DSA",{distance+2},2},
    {"DAD",{distance+3},2},
    {"MMODE",{&monomode},1},{"AMODE",{&analmode},1},
    {"GMODE",{&geomode},1},
    {"GOFFRM",{gammaoff},2},{"GOFFMS",{gammaoff+1},2},
    {"GOFFSA",{gammaoff+2},2},{"GOFFAD",{gammaoff+3},2},
    {"DOFFRM",{deltaoff},2},{"DOFFMS",{deltaoff+1},2},
    {"DOFFSA",{deltaoff+2},2},{"DOFFAD",{deltaoff+3},2},
    {"OFFSETON",{&offseton},1},
    {"XAXIS", {xaxistyp},0},
    {"HC",{dopt},2},{"KC",{dopt+1},2},{"LC",{dopt+2},2},{"EC",{dopt+3},2},
    {"Hi",{hklei},2},{"Ki",{hklei+1},2},{"Li",{hklei+2},2},{"Ei",{hklei+3},2},
    {"H", {hkle},2},{"K", {hkle+1},2},{"L", {hkle+2},2},{"E", {hkle+3},2},
    {"HS",{step},2},{"KS",{step+1},2},{"LS",{step+2},2},{"ES",{step+3},2},
    {"NPTS",{&npts},1},
    {"CSTYPE",{cstyp},0},
    {"H0",{point},2},{"K0",{point+1},2},{"L0",{point+2},2},{"E0",{point+3},2},
    {"H1",{dire1},2},{"K1",{dire1+1},2},{"L1",{dire1+2},2},{"E1",{dire1+3},2},
    {"H2",{dire2},2},{"K2",{dire2+1},2},{"L2",{dire2+2},2},{"E2",{dire2+3},2},
    {"H3",{dire3},2},{"K3",{dire3+1},2},{"L3",{dire3+2},2},{"E3",{dire3+3},2},
    {"HREF1",{ref1},2},{"KREF1",{ref1+1},2},{"LREF1",{ref1+2},2},
    {"HREF2",{ref2},2},{"KREF2",{ref2+1},2},{"LREF2",{ref2+2},2},
    {"HREF3",{ref3},2},{"KREF3",{ref3+1},2},{"LREF3",{ref3+2},2},
    {"reduced_units",{&reduced},1},
    {"AREA", {&sampleArea},2},
    {"THICK",{&sampleThick},2},
    {"FSQ",  {&Fsq},2},
    {"ATTEN",{&absorbCoef},2},
    {"RXFX", {&rxtotalflux},2},
    {"EFXPK",{&epeakflux},2},
    {"MAREA",{&monArea},2},
    {"MEFFK",{&moneff},2},
    {"MREF", {&mref},2},
    {"AREF", {&aref},2},
    {"ATPK", {&calcATpeak},1},
    {"MUOFF", {&muoff}, 2},
    {"IMAX1",{Imax},1},{"IMAX2",{Imax+1},1},{"IMAXZ",{Imax+2},1},
    {"IMAXE",{Imax+3},1},
    {"IMIN1",{Imin},1},{"IMIN2",{Imin+1},1},{"IMINZ",{Imin+2},1},
    {"IMINE",{Imin+3},1},
    {"IREL1",{Irel},2},{"IREL2",{Irel+1},2},{"IRELZ",{Irel+2},2},
    {"IRELE",{Irel+3},2},
    {"IABS1",{Iabs},2},{"IABS2",{Iabs+1},2},{"IABSZ",{Iabs+2},2},
    {"IABSE",{Iabs+3},2},
    {"RPMIN",{&resprobmin},2},
    {(char*)0, {(void*)0}, 0}
  } ;
  
  static ENTRY  outputentries[] = {
    {"HKI",{hkleI},2},{"KKI",{hkleI+1},2},
    {"LKI",{hkleI+2},2},{"EKI",{hkleI+3},2},
    {"HKF",{hkleF},2},{"KKF",{hkleF+1},2},
    {"LKF",{hkleF+2},2},{"EKF",{hkleF+3},2},
    {"TTM",{angls},2},{"OMEGA",{angls+1},2},
    {"TTS",{angls+2},2},{"TTA",{angls+3},2},
    {"TTH0",{&tth0},2},
    {"COTA",{&cota},2},
    {"HPK",{peak},2},{"KPK",{peak+1},2},{"LPK",{peak+2},2},{"EPK",{peak+3},2},
    {"PKMAG",{&peakmag},2},
    {"HFW",{fwhm},2},{"KFW",{fwhm+1},2},{"LFW",{fwhm+2},2},{"EFW",{fwhm+3},2},
    {"FWMAG",{&fwhmmag},2},
    {"MFW",{fwhmang},2},{"WFW",{fwhmang+1},2},
    {"SFW",{fwhmang+2},2},{"AFW",{fwhmang+3},2},
    {"CFW",{&fwhmcot},2},
    {"TOT",{&tot},2},
    {"COR",{&cor},2},
    {"INT",{&relint},2},
    {"CNT",{&maxcountrate},2},
    {"XLN",{&extinctLength},2},
    {"VOL",{&vol},2},
    {"MON",{&mon},2},
    {"RAT",{&hklecountrate},2},
    {"PXYAX",{&prM[2][0][0]},2},{"PXYAY",{&prM[2][0][1]},2},
    {"PXYBX",{&prM[2][1][0]},2},{"PXYBY",{&prM[2][1][1]},2},
    {"PYEAY",{&prM[0][0][0]},2},{"PYEAE",{&prM[0][0][1]},2},
    {"PYEBY",{&prM[0][1][0]},2},{"PYEBE",{&prM[0][1][1]},2},
    {"PEXAE",{&prM[1][0][0]},2},{"PEXAX",{&prM[1][0][1]},2},
    {"PEXBE",{&prM[1][1][0]},2},{"PEXBX",{&prM[1][1][1]},2},
    {"IXYAX",{&scM[2][0][0]},2},{"IXYAY",{&scM[2][0][1]},2},
    {"IXYBX",{&scM[2][1][0]},2},{"IXYBY",{&scM[2][1][1]},2},
    {"IYEAY",{&scM[0][0][0]},2},{"IYEAE",{&scM[0][0][1]},2},
    {"IYEBY",{&scM[0][1][0]},2},{"IYEBE",{&scM[0][1][1]},2},
    {"IEXAE",{&scM[1][0][0]},2},{"IEXAX",{&scM[1][0][1]},2},
    {"IEXBE",{&scM[1][1][0]},2},{"IEXBX",{&scM[1][1][1]},2},
    {"pXYAX",{&prpM[2][0][0]},2},{"pXYAY",{&prpM[2][0][1]},2},
    {"pXYBX",{&prpM[2][1][0]},2},{"pXYBY",{&prpM[2][1][1]},2},
    {"pYEAY",{&prpM[0][0][0]},2},{"pYEAE",{&prpM[0][0][1]},2},
    {"pYEBY",{&prpM[0][1][0]},2},{"pYEBE",{&prpM[0][1][1]},2},
    {"pEXAE",{&prpM[1][0][0]},2},{"pEXAX",{&prpM[1][0][1]},2},
    {"pEXBE",{&prpM[1][1][0]},2},{"pEXBX",{&prpM[1][1][1]},2},
    {"iXYAX",{&scpM[2][0][0]},2},{"iXYAY",{&scpM[2][0][1]},2},
    {"iXYBX",{&scpM[2][1][0]},2},{"iXYBY",{&scpM[2][1][1]},2},
    {"iYEAY",{&scpM[0][0][0]},2},{"iYEAE",{&scpM[0][0][1]},2},
    {"iYEBY",{&scpM[0][1][0]},2},{"iYEBE",{&scpM[0][1][1]},2},
    {"iEXAE",{&scpM[1][0][0]},2},{"iEXAX",{&scpM[1][0][1]},2},
    {"iEXBE",{&scpM[1][1][0]},2},{"iEXBX",{&scpM[1][1][1]},2},
    {"PROJX",{axproj},2},{"PROJY",{axproj+1},2},{"PROJE",{axproj+2},2},
    {"PROJZ",{&vres},2},
    {"PEDAD",{&pedirM[0][0]},2},{"PEDAE",{&pedirM[0][1]},2},
    {"PEDBD",{&pedirM[1][0]},2},{"PEDBE",{&pedirM[1][1]},2},
    {"IEDAD",{&sedirM[0][0]},2},{"IEDAE",{&sedirM[0][1]},2},
    {"IEDBD",{&sedirM[1][0]},2},{"IEDBE",{&sedirM[1][1]},2},
    {"pEDAD",{&pedirpM[0][0]},2},{"pEDAE",{&pedirpM[0][1]},2},
    {"pEDBD",{&pedirpM[1][0]},2},{"pEDBE",{&pedirpM[1][1]},2},
    {"iEDAD",{&sedirpM[0][0]},2},{"iEDAE",{&sedirpM[0][1]},2},
    {"iEDBD",{&sedirpM[1][0]},2},{"iEDBE",{&sedirpM[1][1]},2},
    {"FWHMBX",{fwhmb},2},{"FWHMBY",{fwhmb+1},2},{"FWHMBE",{fwhmb+2},2},
    {"FWHME",{&fwhmeinc},2},{"FWHMZ",{&fwhmz},2},
    {"FWHMPL",{fwhmp},2},{"FWHMPS",{fwhmp+1},2},{"PANGLE",{fwhmp+2},2},
    {"EIG1X",{&eigenvecsM[0][0]},2},{"EIG1Y",{&eigenvecsM[1][0]},2},
    {"EIG1E",{&eigenvecsM[2][0]},2},
    {"EIG2X",{&eigenvecsM[0][1]},2},{"EIG2Y",{&eigenvecsM[1][1]},2},
    {"EIG2E",{&eigenvecsM[2][1]},2},
    {"EIG3X",{&eigenvecsM[0][2]},2},{"EIG3Y",{&eigenvecsM[1][2]},2},
    {"EIG3E",{&eigenvecsM[2][2]},2},
    {"EIG1",{eigenvals},2},{"EIG2",{eigenvals+1},2},{"EIG3",{eigenvals+2},2},
    {"OFFX",{xyeoff},2},{"OFFY",{xyeoff+1},2},{"OFFE",{xyeoff+2},2},
    {"OFFZ",{&qzoff},2},
    {"OFFA",{xyeoffw},2},{"OFFP",{xyeoffw+1},2},
    {"OFFH",{hkloff},2},{"OFFK",{hkloff+1},2},{"OFFL",{hkloff+2},2},
    {"OFFWX",{xyeoffw},2},{"OFFWY",{xyeoffw+1},2},
    {"KIX",{kixy},2},{"KIY",{kixy+1},2},
    {"KIMAG",{&kimag},2},
    {"KFX",{kfxy},2},{"KFY",{kfxy+1},2},
    {"KFMAG",{&kfmag},2},
    {"offX",{&xcnoff},2},{"offY",{&ycnoff},2},{"offZ",{&zcnoff},2},
    {"offE",{&ecnoff},2},{"offKI",{&kioff},2},{"offKF",{&kfoff},2},
    {"RXX",{&rmM[0][0]},2},{"RXY",{&rmM[0][1]},2},{"RXE",{&rmM[0][2]},2},
    {"RYX",{&rmM[1][0]},2},{"RYY",{&rmM[1][1]},2},{"RYE",{&rmM[1][2]},2},
    {"REX",{&rmM[2][0]},2},{"REY",{&rmM[2][1]},2},{"REE",{&rmM[2][2]},2},
    {"RCXX",{&rmcnM[0][0]},2},{"RCXY",{&rmcnM[0][1]},2},
    {"RCXE",{&rmcnM[0][2]},2},
    {"RCYY",{&rmcnM[1][1]},2},{"RCYE",{&rmcnM[1][2]},2},
    {"RCEE",{&rmcnM[2][2]},2},
    {"RIXX",{&rminv[0][0]},2},{"RIXY",{&rminv[0][1]},2},
    {"RIXE",{&rminv[0][2]},2},
    {"RIYY",{&rminv[1][1]},2},{"RIYE",{&rminv[1][2]},2},
    {"RIEE",{&rminv[2][2]},2},
    {"DET",{&detrm},2},
    {"GCXX",{&gc[0][0]},2},{"GCXY",{&gc[0][1]},2},{"GCXE",{&gc[0][2]},2},
    {"GCYY",{&gc[1][1]},2},{"GCYE",{&gc[1][2]},2},{"GCEE",{&gc[2][2]},2},
    {"GLX",{gl},2},{"GLY",{gl+1},2},{"GLE",{gl+2},2},
    {"DMDX",{&dfh[0][0]},2},{"DMDY",{&dfh[0][1]},2},{"DMDE",{&dfh[0][2]},2},
    {"DMFX",{&dfh[1][0]},2},{"DMFY",{&dfh[1][1]},2},{"DMFE",{&dfh[1][2]},2},
    {"DMHX",{&dfh[2][0]},2},{"DMHY",{&dfh[2][1]},2},{"DMHE",{&dfh[2][2]},2},
    {"BRMX",{brm},2},{"BRMY",{brm+1},2},{"BRME",{brm+2},2},
    {"C2XX",{&c2bar[0][0]},2},{"C2XY",{&c2bar[0][1]},2},
    {"C2XE",{&c2bar[0][2]},2},
    {"C2YY",{&c2bar[1][1]},2},{"C2YE",{&c2bar[1][2]},2},
    {"C2EE",{&c2bar[2][2]},2},
    {"BB1",{bbar},2},{"BB2",{bbar+1},2},{"BB3",{bbar+2},2},
    {"A1",{avec},2},{"A2",{avec+1},2},{"A3",{avec+2},2},
    {"B1",{bvec},2},{"B2",{bvec+1},2},{"B3",{bvec+2},2},
    {"C1",{cvec},2},{"C2",{cvec+1},2},{"C3",{cvec+2},2},
    {"Z1",{zuni},2},{"Z2",{zuni+1},2},{"Z3",{zuni+2},2},
    {"AST1",{ast},2},{"AST2",{ast+1},2},{"AST3",{ast+2},2},
    {"BST1",{bst},2},{"BST2",{bst+1},2},{"BST3",{bst+2},2},
    {"CST1",{cst},2},{"CST2",{cst+1},2},{"CST3",{cst+2},2},
    {"AAST",{dpstar},2},{"ABST",{dpstar+1},2},{"ACST",{dpstar+2},2},
    {"BBST",{dpstar+3},2},{"BCST",{dpstar+4},2},{"CCST",{dpstar+5},2},
    {"RLAT1",{rlat},2},{"RLAT2",{rlat+1},2},{"RLAT3",{rlat+2},2},
    {"RANG1",{rang},2},{"RANG2",{rang+1},2},{"RANG3",{rang+2},2},
    {"AUN1",{auni},2},{"AUN2",{auni+1},2},{"AUN3",{auni+2},2},
    {"AUNX1",{aunix},2},{"AUNX2",{aunix+1},2},{"AUNX3",{aunix+2},2},
    {"AHKLX1",{ahklx},2},{"AHKLX2",{ahklx+1},2},{"AHKLX3",{ahklx+2},2},
    {"AMAG",{&amag},2},
    {"APU1",{apun},2},{"APU2",{apun+1},2},{"APU3",{apun+2},2},
    {"APUX1",{apunx},2},{"APUX2",{apunx+1},2},{"APUX3",{apunx+2},2},
    {"BHKLX1",{bhklx},2},{"BHKLX2",{bhklx+1},2},{"BHKLX3",{bhklx+2},2},
    {"BMAG",{&bmag},2},
    {"ABDOT",{&abdot},2},
    {"APBDOT",{&apbdot},2},
    {"APB1",{apub},2},{"APB2",{apub+1},2},{"APB3",{apub+2},2},
    {"APBH",{phkln},2},{"APBK",{phkln+1},2},{"APBL",{phkln+2},2},
    {"DIRM",{spin},1},{"DIRS",{spin+1},1},{"DIRA",{spin+2},1},
    {"X0",{q0},2},{"Y0",{q0+1},2},
    {"REF1X",{ref1xy},2},{"REF1Y",{ref1xy+1},2},
    {"REF1MAG",{&ref1mag},2},
    {"REF2X",{ref2xy},2},{"REF2Y",{ref2xy+1},2},
    {"REF2MAG",{&ref2mag},2},
    {"REF3X",{ref3xy},2},{"REF3Y",{ref3xy+1},2},
    {"REF3MAG",{&ref3mag},2},
    {"XCSPT",{cspoint},2},{"YCSPT",{cspoint+1},2},{"ECSPT",{cspoint+3},2},
    {"XCSD1",{csdire1},2},{"YCSD1",{csdire1+1},2},{"ECSD1",{csdire1+3},2},
    {"XCSD2",{csdire2},2},{"YCSD2",{csdire2+1},2},{"ECSD2",{csdire2+3},2},
    {"XCSU1",{cs1diru},2},{"YCSU1",{cs1diru+1},2},
    {"QMAG",{&hklmag},2},{"iMAG",{&ihklmag},2},
    {"SMAG",{&stepmag},2},
    {"PTMAG",{&pointmag},2},
    {"D1MAG",{&diremag[0]},2},{"D2MAG",{&diremag[1]},2},
    {"D3MAG",{&diremag[2]},2},
    {"HQPU",{hklqpun},2},{"KQPU",{hklqpun+1},2},{"LQPU",{hklqpun+2},2},
    {"QPMAG",{&qpmag},2},
    {"EFIX",{&efix},2},
    {"NEIG",{&nfree},1},
    {(char*)0, {(void*)0}, 0}
  } ;
