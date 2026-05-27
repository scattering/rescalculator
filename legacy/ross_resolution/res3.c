/*
 * 2 and 3 axis resolution and convolution functions
 * for the lsq3DAVEinterface
 * this is basically resRes.c from the tclres package
 * with the tcl stuff removed.
 *
 *    This module implements "neutron spectrometer resolution"
 *    and convolution preparation
 *
 * since this will become a load module
 * make sure local memory allocation is deallocated
 * provide an entrypoint for memory cleanup = Res3_shutdown()
 */


/* global static res in and out variables */

#include "res3int.h"




static void startup() ;

int Res3_Init()
{
  /* call startup to initialize all the dictionaries etc */
  startup() ;
}

static int setErrmsg(char *msg)
{
  int nc ;
  if( ! msg ) return 0 ;
  if( (nc = strlen(msg)) >= errmsgAlloc ) {
    errmsg = (char*)realloc(errmsg, (nc+1)*sizeof(char)) ;
    if( ! errmsg ) return 0 ;
    errmsgAlloc = nc+1 ;
  }
  strcpy(errmsg, msg) ;
  return 1 ;
}
void Res3_clearErrmsg()
{
  if( errmsg ) errmsg[0] = '\0' ;
}
char *Res3_getErrmsg()
{
  return errmsg ;
}

static int recipCalc(double *lat, double *ang,
		     double *av, double *bv, double *cv,
		     double *as, double *bs, double *cs, double *cvol) ;


static void dores3calc() ;
static void dores3rm() ;
static int dores3rmab() ;


static int Res3_CalcOp()		       
{
  return dores3rmab() ;
}

int Res3_GetInput(char *name, char **data)
{
  /*
    get the data for input item name
  */
  ENTRY *entryptr ;
  ENTRY entryA ;

  entryA.type = -1 ;  // flags fetch entry data
  entryA.name = name ;
  if( ! (entryptr = dictionary( globalinputdict, &entryA, 0 )) ) return 0 ;
  if ( entryptr->type == 2 && fabs(*(entryptr->data.d)) < ISZERO )
    *(entryptr->data.d) = 0. ;
  if( ! entrydatatostringNL( entryptr, data ) ) return 0 ;
  return 1 ;
}
int Res3_GetNextInput(char **name, char **data, int init)
{
  /*
    get the next name,data pair from input dictionary
  */
  return dictgetnextentry(globalinputdict, name, data, init) ;
}
int Res3_GetNextOutput(char **name, char **data, int init)
{
  /*
    get the next name,data pair from output dictionary
  */
  return dictgetnextentry(globaloutputdict, name, data, init) ;
}

int Res3_GetOutput(char *name, char **data)
{
  /*
    get the data for output item name
  */
  ENTRY *entryptr ;
  ENTRY entryA ;

  entryA.type = -1 ;  // flags fetch entry data
  entryA.name = name ;
  if( ! (entryptr = dictionary( globaloutputdict, &entryA, 0 )) ) return 0 ;
  if ( entryptr->type == 2 && fabs(*(entryptr->data.d)) < ISZERO )
    *(entryptr->data.d) = 0. ;
  if( ! entrydatatostringNL( entryptr, data ) ) return 0 ;
  return 1 ;
}



int Res3_SetOp(int argc, char **argv)
{
  /*
   * sets  symbol1  value1  symbol2  value2 ...
   */

  int i ;

  if ( argc < 1 || argv == NULL ) return 0 ;

  for( i=1 ; i<argc ; i += 2 )
    {
      stringtodata( globalinputdict, argv[i-1], argv[i] ) ;
    }
  return 1 ;
}



/* dictionary precision */
void Res3_setDblprec( int prec )
{
  setDblprec(prec) ;
}




static int Res3_WriteOp(int typ, char *fil)
{
  /*
   * typ=0 inputvalues  typ=1 outputvalues  typ=2 both
   */


  ENTRY *entryptr ;
  char buf[256] ;

  FILE *fpt ;

  if ( ! fil ) return 0 ;

  if( (fpt = fopen(fil, "w")) == (FILE*)NULL ) {
    return 0 ;
  }

  if( ! typ%2 ) {
    entryptr = globalinput ;
    while( entryptr->name != (char*)0 )
      {
	if ( entryptr->type == 2 )
	  if( fabs(*(entryptr->data.d)) < ISZERO )
	    *(entryptr->data.d) = 0. ;
	entrydatatostringN( entryptr, buf ) ;
	fprintf(fpt, "%s %s\n", entryptr->name, buf) ;
	++entryptr ;
      }
  }
  if( typ ) {
    entryptr = globaloutput ;
    while( entryptr->name != (char*)0 )
      {
	if ( entryptr->type == 2 )
	  if( fabs(*(entryptr->data.d)) < ISZERO )
	    *(entryptr->data.d) = 0. ;
	entrydatatostringN( entryptr, buf ) ;
	fprintf(fpt, "%s %s\n", entryptr->name, buf) ;
	++entryptr ;
      }
  }
  fclose(fpt) ;
  return 1 ;
}



int Res3_checkConvNpts(int np, Conv3_conv3 *conv3pt)
{
  int i ;
  Conv3_resinfo *respt ;
  if ( np > conv3pt->nalloc ) {
    conv3pt->convpts =
      (Conv3_convpt *)realloc(conv3pt->convpts, np*sizeof(Conv3_convpt)) ;
    if( ! conv3pt->convpts ) return 0 ;
    for( i=conv3pt->nalloc ; i<np ; i++ ) {
      conv3pt->convpts[i].nres = 0 ;
      conv3pt->convpts[i].nalloc = 0 ;
      conv3pt->convpts[i].respts = NULL ;
    }
    conv3pt->y = (double *)realloc(conv3pt->y, np*sizeof(double)) ;
    if( ! conv3pt->y ) return 0 ;
    conv3pt->ye = (double *)realloc(conv3pt->ye, np*sizeof(double)) ;
    if( ! conv3pt->ye ) return 0 ;
    conv3pt->h = (double *)realloc(conv3pt->h, np*sizeof(double)) ;
    if( ! conv3pt->h ) return 0 ;
    conv3pt->k = (double *)realloc(conv3pt->k, np*sizeof(double)) ;
    if( ! conv3pt->k ) return 0 ;
    conv3pt->l = (double *)realloc(conv3pt->l, np*sizeof(double)) ;
    if( ! conv3pt->l ) return 0 ;
    conv3pt->e = (double *)realloc(conv3pt->e, np*sizeof(double)) ;
    if( ! conv3pt->e ) return 0 ;
    conv3pt->Ei = (double *)realloc(conv3pt->Ei, np*sizeof(double)) ;
    if( ! conv3pt->Ei ) return 0 ;
    conv3pt->Ef = (double *)realloc(conv3pt->Ef, np*sizeof(double)) ;
    if( ! conv3pt->Ef ) return 0 ;
    conv3pt->fs = (double *)realloc(conv3pt->fs, np*sizeof(double)) ;
    if( ! conv3pt->fs ) return 0 ;
    conv3pt->omega = (double *)realloc(conv3pt->omega, np*sizeof(double)) ;
    if( ! conv3pt->omega ) return 0 ;
    conv3pt->chiA = (double *)realloc(conv3pt->chiA, np*sizeof(double)) ;
    if( ! conv3pt->chiA ) return 0 ;
    conv3pt->chiB = (double *)realloc(conv3pt->chiB, np*sizeof(double)) ;
    if( ! conv3pt->chiB ) return 0 ;
    conv3pt->angleA = (double *)realloc(conv3pt->angleA, np*sizeof(double)) ;
    if( ! conv3pt->angleA ) return 0 ;
    conv3pt->T = (double *)realloc(conv3pt->T, np*sizeof(double)) ;
    if( ! conv3pt->T ) return 0 ;
    conv3pt->H = (double *)realloc(conv3pt->H, np*sizeof(double)) ;
    if( ! conv3pt->H ) return 0 ;
    conv3pt->Pz = (double *)realloc(conv3pt->Pz, np*sizeof(double)) ;
    if( ! conv3pt->Pz ) return 0 ;
    conv3pt->Px = (double *)realloc(conv3pt->Px, np*sizeof(double)) ;
    if( ! conv3pt->Px ) return 0 ;
    conv3pt->Py = (double *)realloc(conv3pt->Py, np*sizeof(double)) ;
    if( ! conv3pt->Py ) return 0 ;
    for( i=conv3pt->nalloc ; i<np ; i++ ) {
      conv3pt->y[i] = 0. ;
      conv3pt->ye[i] = 0. ;
      conv3pt->h[i] = 0. ;
      conv3pt->k[i] = 0. ;
      conv3pt->l[i] = 0. ;
      conv3pt->e[i] = 0. ;
      conv3pt->Ei[i] = 0. ;
      conv3pt->Ef[i] = 0. ;
      conv3pt->fs[i] = 0. ;
      conv3pt->omega[i] = 0. ;
      conv3pt->chiA[i] = 0. ;
      conv3pt->chiB[i] = 0. ;
      conv3pt->angleA[i] = 0. ;
      conv3pt->T[i] = 0. ;
      conv3pt->H[i] = 0. ;
      conv3pt->Pz[i] = 0. ;
      conv3pt->Px[i] = 0. ;
      conv3pt->Py[i] = 0. ;
    }
    conv3pt->nalloc = np ;
  }
  if( np > conv3pt->ndata ) conv3pt->ndata = np ;
  respt = conv3pt->res ;
  if( respt != NULL ) respt->npts = conv3pt->ndata ;
  return 1 ;
}

static int
Conv3_NptsOp(Conv3_conv3 *conv3pt, int np)
{

  if( conv3pt == NULL ) return 0 ;
  if ( np <= 0 ) return conv3pt->ndata ;
  if( ! Res3_checkConvNpts(np, conv3pt) ) return 0 ;
  conv3pt->ndata = np ;
  return np ;
}

Conv3_resinfo *Res3_NewResinfo() ;
static int Copyconv3(Conv3_conv3 *destPtr, Conv3_conv3 *srcPtr) ;


Conv3_conv3 *Res3_Conv3MeshPrep(Conv3_conv3 *conv3pt,
				Conv3_conv3 *interconv, int mesh)
{
  /*
   * prepare mesh data before calc
   */

  int i, j, ix ;
  int mesh1, meshpts ;
  double jm, dmesh1 ;
  double dh, dk, dl, de, dEi, dEf, dT, dH, dPx, dPy, dPz ;
  double domega, dchiA, dchiB, dangleA ;

  if( ! conv3pt || ! conv3pt->res ) return NULL ;
  Res3_clearErrmsg() ;
  if ( conv3pt->ndata < 1 ) {
    setErrmsg("ERROR: convMeshPrep source conv has no data!") ;
    return NULL ;
  }
  if( ! interconv ) interconv = Res3_NewConv3s() ;
  if( ! interconv ) return NULL ;

  if( mesh < 0 ) mesh = 0 ;
  mesh1 = mesh + 1 ;
  dmesh1 = mesh1 ;
  Copyconv3(interconv, conv3pt) ;
  meshpts = (conv3pt->ndata - 1)*mesh1 + 1 ;
  if( ! Res3_checkConvNpts(meshpts, interconv) ) {
    setErrmsg(res3memerr) ;
    return NULL ;
  }

  ix = 0 ;
  for( i=0 ; i<conv3pt->ndata - 1 ; i++ ) {
    dh = conv3pt->h[i+1] - conv3pt->h[i] ;
    dk = conv3pt->k[i+1] - conv3pt->k[i] ;
    dl = conv3pt->l[i+1] - conv3pt->l[i] ;
    de = conv3pt->e[i+1] - conv3pt->e[i] ;
    dEi = conv3pt->Ei[i+1] - conv3pt->Ei[i] ;
    dEf = conv3pt->Ef[i+1] - conv3pt->Ef[i] ;
    domega = conv3pt->omega[i+1] - conv3pt->omega[i] ;
    dchiA = conv3pt->chiA[i+1] - conv3pt->chiA[i] ;
    dchiB = conv3pt->chiB[i+1] - conv3pt->chiB[i] ;
    dangleA = conv3pt->angleA[i+1] - conv3pt->angleA[i] ;
    dT = conv3pt->T[i+1] - conv3pt->T[i] ;
    dH = conv3pt->H[i+1] - conv3pt->H[i] ;
    dPx = conv3pt->Px[i+1] - conv3pt->Px[i] ;
    dPy = conv3pt->Py[i+1] - conv3pt->Py[i] ;
    dPz = conv3pt->Pz[i+1] - conv3pt->Pz[i] ;
    for( j=0 ; j<mesh1 ; j++ ) {
      jm = j/dmesh1 ;
      interconv->y[ix] = conv3pt->y[i] ;
      interconv->ye[ix] = conv3pt->ye[i] ;
      interconv->fs[ix] = conv3pt->fs[i] ;
      interconv->h[ix] = conv3pt->h[i] + jm*dh ;
      interconv->k[ix] = conv3pt->k[i] + jm*dk ;
      interconv->l[ix] = conv3pt->l[i] + jm*dl ;
      interconv->e[ix] = conv3pt->e[i] + jm*de ;
      interconv->Ei[ix] = conv3pt->Ei[i] + jm*dEi ;
      interconv->Ef[ix] = conv3pt->Ef[i] + jm*dEf ;
      interconv->T[ix] = conv3pt->T[i] + jm*dT ;
      interconv->H[ix] = conv3pt->H[i] + jm*dH ;
      interconv->Px[ix] = conv3pt->Px[i] + jm*dPx ;
      interconv->Py[ix] = conv3pt->Py[i] + jm*dPy ;
      interconv->Pz[ix] = conv3pt->Pz[i] + jm*dPz ;

      interconv->omega[ix] = conv3pt->omega[i] + jm*domega ;
      interconv->chiA[ix] = conv3pt->chiA[i] + jm*dchiA ;
      interconv->chiB[ix] = conv3pt->chiB[i] + jm*dchiB ;
      interconv->angleA[ix] = conv3pt->angleA[i] + jm*dangleA ;

      ix++ ;
    }
  }
  i = conv3pt->ndata - 1 ;
  interconv->y[ix] = conv3pt->y[i] ;
  interconv->ye[ix] = conv3pt->ye[i] ;
  interconv->fs[ix] = conv3pt->fs[i] ;
  interconv->h[ix] = conv3pt->h[i] ;
  interconv->k[ix] = conv3pt->k[i] ;
  interconv->l[ix] = conv3pt->l[i] ;
  interconv->e[ix] = conv3pt->e[i] ;
  interconv->Ei[ix] = conv3pt->Ei[i] ;
  interconv->Ef[ix] = conv3pt->Ef[i] ;
  interconv->T[ix] = conv3pt->T[i] ;
  interconv->H[ix] = conv3pt->H[i] ;
  interconv->Px[ix] = conv3pt->Px[i] ;
  interconv->Py[ix] = conv3pt->Py[i] ;
  interconv->Pz[ix] = conv3pt->Pz[i] ;
  
  interconv->omega[ix] = conv3pt->omega[i] ;
  interconv->chiA[ix] = conv3pt->chiA[i] ;
  interconv->chiB[ix] = conv3pt->chiB[i] ;
  interconv->angleA[ix] = conv3pt->angleA[i] ;
  interconv->ndata = meshpts ;
  if( ! Res3_Conv3Prep(interconv) ) return NULL ;
  return interconv ;
}



static void symeig2(double pxx, double pyy, double pxy,
		    double *eig1, double *eig1x, double *eig1y,
		    double *eig2, double *eig2x, double *eig2y)
{
  double thet, st, ct, st2, ct2, sc2 ;

  thet = 0.5*atan2(2.*pxy, pyy-pxx) ;
  st = sin(thet) ;
  ct = cos(thet) ;
  *eig1x = ct ;
  *eig1y = -st ;
  *eig2x = st ;
  *eig2y = ct ;
  st2 = st*st ;
  ct2 = ct*ct ;
  sc2 = 2.*st*ct ;

  *eig1 = ct2*pxx - sc2*pxy + st2*pyy ;
  *eig2 = st2*pxx + sc2*pxy + ct2*pyy ;
}

static int resprojcalc( double *r, Conv3_respt *ires )
{
  /* NO globals required */

  double eig1, eig1x, eig1y, eig2, eig2x, eig2y ;
  double rxe, rye ;
  double detxy, x0el, y0el ;
  double pxx, pyy, pxy ;

  if( r[5] <= 0. ) r[5] = 1.e8 ;

  /*
    For elastic cross sections
    Also need to get intersection with E=0 plane which will
    move Q0 center and change the value of ratmin1 below
  */

  detxy = r[0]*r[3] - r[1]*r[1] ;
  if( detxy <= 0. ) { return 0 ; }
  x0el = ires->E0*(r[1]*r[4] - r[3]*r[2])/detxy ;
  y0el = ires->E0*(r[1]*r[2] - r[0]*r[4])/detxy ;
  ires->Q0el[0] = ires->Q0[0] - x0el ;
  ires->Q0el[1] = ires->Q0[1] - y0el ;
  ires->Q0el[2] = ires->Q0[2] ;
  //ires->Q0[0] -= x0el ;
  //ires->Q0[1] -= y0el ;
  /* also compute correction ro resprob */
  ires->rprob = r[5]*ires->E0*ires->E0 +
    r[0]*x0el*x0el + 2.*r[1]*x0el*y0el + r[3]*y0el*y0el +
    2.*ires->E0*(r[2]*x0el + r[4]*y0el) ;
  pxx = r[0] ;
  pyy = r[3] ;
  pxy = r[1] ;

  symeig2(pxx, pyy, pxy, &eig1, &eig1x, &eig1y, &eig2, &eig2x, &eig2y) ;

  if( eig1 <= 0. || eig2 <= 0. ) { return 0 ; }
  

  ires->ieigA[0] = eig1x ;
  ires->ieigA[1] = eig1y ;
  ires->ieigB[0] = eig2x ;
  ires->ieigB[1] = eig2y ;
  ires->ieiga = eig1 ;
  ires->ieigb = eig2 ;
 
 
  /* now for inelastic use projection onto E plane */
  
  rxe = -r[2]/r[5] ;
  rye = -r[4]/r[5] ;
  ires->eoptXY[0] = rxe ;
  ires->eoptXY[1] = rye ;
  pxx = r[0] + rxe*r[2] ;
  pyy = r[3] + rye*r[4] ;
  pxy = r[1] + rxe*r[4] ;

  symeig2( pxx, pyy, pxy, &eig1, &eig1x, &eig1y, &eig2, &eig2x, &eig2y ) ;

  if( eig1 <= 0. || eig2 <= 0. ) { return 0 ; }
  
  ires->peigA[0] = eig1x ;
  ires->peigA[1] = eig1y ;
  ires->peigB[0] = eig2x ;
  ires->peigB[1] = eig2y ;
  ires->peiga = eig1 ;
  ires->peigb = eig2 ;  

  if( r[6] <= 0. ) r[6] = 1.e9 ;
  ires->nz = 1./sqrt(r[6]) ;
  return 1 ;
}




static void copy33( M33 src, M33 dest )
{
  int i, j ;
  for( i=0 ; i<3 ; i++ ) {
    for( j=0 ; j<3 ; j++ ) {
      dest[i][j] = src[i][j] ;
    }
  }
}
static void transform(M33 T, double *in, double *out)
{
  out[0] = T[0][0]*in[0] + T[0][1]*in[1] + T[0][2]*in[2] ;
  out[1] = T[1][0]*in[0] + T[1][1]*in[1] + T[1][2]*in[2] ;
  out[2] = T[2][0]*in[0] + T[2][1]*in[1] + T[2][2]*in[2] ;
}

Conv3_resinfo *Res3_NewResinfo()
{
  Conv3_resinfo *respt ;
  static Conv3_resinfo defltResinfo =
    {
      {4., 4., 4.}, {90.,90.,90.},
      {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.},
      {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.},
      {1., 0., 0.}, {0., 0., 1.},
      {1., 0., 0.}, {0., 0., 1.},
      {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.},
      {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.},
      1., 1.,
      { {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.} },
      { {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.} },
      { {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.} },
      { {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.} },
      { {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.} },
      {40., 40., 40., 40.}, {180., 180., 180., 180.},
      {1.,   2.,  2.,  4.}, {180., 180., 180., 180.},
      {40., 40., 40.}, {40., 40., 40.},
      {3.35416, 3.35416},
      {1., 1., 1., 1.},
      14.7, 0., 2.351,
      0., 0., 0.,
      0, 0, 0, 0, 1, 0, 0, {1, -1, 1},
      "rlr", "ahkl", NULL, NULL, NULL,
      0.04, 35., 24., 1.e-5, 1., 1.,
      { 8, 8, 8, 8 },
      { 0, 0, 0, 0 },
      {1.e-2, 1.e-2, 1.e-2, 1.e-2},
      {1.e-6, 1.e-6, 1.e-6, 1.e-6},
      1.e-4,
      {300., 0., 0., 0., 0.}
    } ;

  respt = (Conv3_resinfo *)calloc(1, sizeof(Conv3_resinfo)) ;
  if( ! respt ) return 0 ;
  *respt = defltResinfo ;
  return respt ;
}

static mxFun mxLookup(char *id)
{
  int i ;
  if( ! id ) return NULL ;
  for( i=0 ; i<nmxfunctions ; i++ ) {
    if( strstr(mxfunctions[i].uid, id) ) return mxfunctions[i].f ;
  }
  return NULL ;
}

int Res3_Conv3Prep(Conv3_conv3 *conv3pt)
{
  /*
   *
   *  Use conv3 data and resinfo to calc all conv points
   *  The xtalprogs must return number of xtal elements with no arg call
   *  and with index arg calc the geometry for that indexed sub xtal
   *  using things like  res3 set OMOFFM omegaoff
   */

  int i, j, k, m, code ;
  int nXtal, iXtal ;
  int multiMono, multiAnal, nMono, nAnal ;
  double angs[6] ;
  char jbuf[16], kbuf[16] ;
  char errbuf[128] ;
  static double Dn = 2.072141789 ;

  Conv3_convpt *convpt ;
  Conv3_respt *ires ;
  Conv3_resinfo *respt ;
  static Conv3_resinfo *resave = NULL ;

  int (*mxMono)() ;
  int (*mxAnal)() ;

  if( ! conv3pt ) return 0 ;
  convpt = conv3pt->convpts ;
  respt = conv3pt->res ;
  if( ! convpt || ! respt ) return 0 ;
  Res3_clearErrmsg() ;

  if( ! resave ) {
    if( ! (resave = Res3_NewResinfo()) ) {
      setErrmsg(res3memerr) ;
      return 0 ;
    }
  }


  /*
    save current global resinfo
    and replace with that from this resinfo
  */

  Res3_loadResinfo(resave) ;
  Res3_unloadResinfo(respt) ;

  multiMono = 0 ;
  multiAnal = 0 ;
  nMono = 1 ;
  nAnal = 1 ;

  if( respt->monoprog != NULL ) {
    if( (mxMono = mxLookup(respt->monoprog)) == NULL ) {
      if( nMono = mxMono(0, 0.) ) {
	multiMono = 1 ;
      }
    }
  }
  if( respt->analprog != NULL ) {
    if( (mxAnal = mxLookup(respt->analprog)) == NULL ) {
      if( nAnal = mxAnal(0, 0.) ) {
	multiAnal = 1 ;
      }
    }
  }
  nXtal = nMono*nAnal ;

  /* make sure any pre-existing points start with zero res xtal pairs */
  for( i=0 ; i<conv3pt->ndata ; i++ ) {
    if( nXtal > convpt[i].nalloc ) {
      convpt[i].respts =
	(Conv3_respt *)realloc(convpt[i].respts, nXtal*sizeof(Conv3_respt)) ;
      if( ! convpt[i].respts ) {
	setErrmsg(res3memerr) ;
	return 0 ;
      }
      convpt[i].nalloc = nXtal ;
    }
    convpt[i].nres = nXtal ;

    /* set the global resinfo */
    chiA = conv3pt->chiA[i] ;
    chiB = conv3pt->chiB[i] ;
    angleA = conv3pt->angleA[i] ;
    dopt[0] = conv3pt->h[i] ;
    dopt[1] = conv3pt->k[i] ;
    dopt[2] = conv3pt->l[i] ;
    dopt[3] = conv3pt->e[i] ;
    Ki = sqrt(conv3pt->Ei[i]/Dn) ;
    Kf = sqrt(conv3pt->Ef[i]/Dn) ;

    /* set up the coord system */
    if( ! set_recip() ) {
      Res3_unloadResinfo(resave) ;
      sprintf(errbuf, "Res3_ConvPrep failed set_recip at datapoint %d", i) ;
      setErrmsg(errbuf) ;
      return 0 ;
    }
    if( ! sampangles(dopt, ast, bst, cst, aunix, apub,
		     conv3pt->Ei[i], conv3pt->Ef[i], angs) ) {
      Res3_unloadResinfo(resave) ;
      sprintf(errbuf, "Res3_ConvPrep failed sampangles at datapoint %d", i) ;
      setErrmsg(errbuf) ;
      return 0 ;
    }

    convpt[i].omega = angs[3] ;
    convpt[i].Q12toQantibeam[0] = sin(convpt[i].omega) ;
    convpt[i].Q12toQperpbeam[1] = sin(convpt[i].omega) ;
    if( scattsideright ) {
      convpt[i].Q12toQantibeam[1] = -cos(convpt[i].omega) ;
      convpt[i].Q12toQperpbeam[0] = cos(convpt[i].omega) ;      
    } else {
      convpt[i].Q12toQantibeam[1] = cos(convpt[i].omega) ;
      convpt[i].Q12toQperpbeam[0] = -cos(convpt[i].omega) ;
    }

    copy33( MX, convpt[i].QtoHKL ) ;
    copy33( MI, convpt[i].HKLtoQ ) ;
    copy33( L, convpt[i].L ) ;

    /* also copy the Irel and Iabs from resinfo into each convpt */
    for( j=0 ; j<4 ; j++ ) {
      convpt[i].Imax[j] = respt->Imax[j] ;
      convpt[i].Imin[j] = respt->Imin[j] ;
      convpt[i].Irel[j] = respt->Irel[j] ;
      convpt[i].Iabs[j] = respt->Iabs[j] ;
    }

    if( ! dores3rmab() ) {
      Res3_unloadResinfo(resave) ;
      return 0 ;
    }
    //code = Tcl_Eval(interp, "res3 calc") ;


    /*
      multixtal procs
      best way is to set all geom stuff and turn on offsets
      then calcresm

      or calc  hkleff kie kfe
      and geo collimation stuff
      then use these to do calcresm with offsets off
      this doesnt get it quite right (e.g. coll + geocoll effects and r0)
    */

    iXtal = 0 ;
    for( j=1 ; j<=nMono ; j++ ) {
      if( multiMono ) {
	if( ! mxMono(j, Ki) ) {
	  Res3_unloadResinfo(resave) ;
	  setErrmsg("ERROR: failed mxMono") ;
	  return 0 ;
	}
      }
      for( k=1 ; k<=nAnal ; k++ ) {
	if( multiAnal ) {
	  if( ! mxAnal(k, Kf) ) {
	    setErrmsg("ERROR: failed mxAnal") ;
	    Res3_unloadResinfo(resave) ;
	    return 0 ;
	  }
	}
	if( ! dores3rmab() ) {
	  Res3_unloadResinfo(resave) ;
	  return 0 ;
	}
	/* now copy the calc res stuff to ires */
	ires = convpt[i].respts + iXtal ;
	ires->i = i ;
	for( m=0 ; m<3 ; m++ ) ires->hkl[m] = hkleff[m] ;
	ires->E0 = hkleff[3] ;
	ires->EI = Dn*kie*kie ;
	ires->EF = Dn*kfe*kfe ;
	ires->rm[0] = rm[0][0] ;
	ires->rm[1] = rm[0][1] ;
	ires->rm[2] = rm[0][2] ;
	ires->rm[3] = rm[1][1] ;
	ires->rm[4] = rm[1][2] ;
	ires->rm[5] = rm[2][2] ;
	ires->rm[6] = rmv ;
	ires->r0 = r0 ;
	ires->vol = vol ;

	transform( MI, ires->hkl, ires->Q0 ) ;
	if( ! resprojcalc( ires->rm, ires ) ) {
	  Res3_unloadResinfo(resave) ;
	  return 0 ;
	}
	iXtal++ ;
      }
    }
  }
  Res3_unloadResinfo(resave) ;
  return 1 ;
}

static int hexAddrOK(char *buf)
{
  int nc ;
  char *cp ;

  nc = strlen(buf) ;
  if( nc < 4 ) return 0 ;
  cp = buf ;
  if( buf[0] == '0' && (buf[1] == 'x' || buf[1] == 'X') ) cp = buf + 2 ;
  if( strlen(cp) < 4 ) return 0 ;
  while( *cp != '\0' ) {
    if( *cp < 48 ) return 0 ;
    if( *cp > 102) return 0 ;
    if( *cp > 57 && *cp < 65 ) return 0 ;
    if( *cp > 70 && *cp < 97 ) return 0 ;
    cp++ ;
  }
  return 1 ;
}




/* manage storage for global storage arrays and fill from string */
static int getdbls( char *s )	  
{
  int i ;
  char *cp ;
  double valu ;
  if( (cp = strtok(s, " ,")) == NULL ) return (0) ;
  i = 0 ;
  while( cp != NULL && sscanf(cp, "%lf", &valu) == 1 )
    {
      if( i >= ndblalloc )
	{
	  ndblalloc += incralloc ;
	  if( (dbls = realloc(dbls, ndblalloc*sizeof(double))) == NULL )
	    {
	      ndblalloc = 0 ;
	      return (0) ;
	    }
	}
      dbls[i++] = valu ;
      cp = strtok(NULL," ,") ;
    }
  return (i) ;
}

static int set_recipR(Conv3_resinfo *rp) ;

/*
  readConvFile gets replaced by the reading functions
  in lsq3fit for DAVE interface
*/

Conv3_conv3 *Res3_NewConv3s()
{
  static Conv3_conv3 newConv3s = 
    {
      0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
      NULL, NULL, NULL, NULL, NULL,
      NULL, NULL, NULL, NULL, NULL,
      NULL, NULL,
      NULL
    } ;

  Conv3_conv3 *lPtr ;

  lPtr = (Conv3_conv3 *)calloc(1, sizeof(Conv3_conv3)) ;
  if( ! lPtr ) return 0 ;
  *lPtr = newConv3s ;

  lPtr->res = Res3_NewResinfo() ;

  return lPtr ;
}


void Res3_FreeConv3s(Conv3_conv3 *lPtr)
{
  int i ;
  Conv3_convpt *pPtr ;
  if( ! lPtr ) return ;
  pPtr = lPtr->convpts ;
  if ( pPtr != NULL ) {
    for( i=0 ; i<lPtr->nalloc ; i++ ) {
      if ( pPtr[i].respts != NULL ) free(pPtr[i].respts) ;
    }
    free(pPtr) ;
  }

  if ( lPtr->y != NULL ) free(lPtr->y) ;
  if ( lPtr->ye != NULL ) free(lPtr->ye) ;
  if ( lPtr->h != NULL ) free(lPtr->h) ;
  if ( lPtr->k != NULL ) free(lPtr->k) ;
  if ( lPtr->l != NULL ) free(lPtr->l) ;
  if ( lPtr->e != NULL ) free(lPtr->e) ;
  if ( lPtr->Ei != NULL ) free(lPtr->Ei) ;
  if ( lPtr->Ef != NULL ) free(lPtr->Ef) ;
  if ( lPtr->fs != NULL ) free(lPtr->fs) ;
  if ( lPtr->omega != NULL ) free(lPtr->omega) ;
  if ( lPtr->chiA != NULL ) free(lPtr->chiA) ;
  if ( lPtr->chiB != NULL ) free(lPtr->chiB) ;
  if ( lPtr->angleA != NULL ) free(lPtr->angleA) ;

  if ( lPtr->T != NULL ) free(lPtr->T) ;
  if ( lPtr->H != NULL ) free(lPtr->H) ;
  if ( lPtr->Pz != NULL ) free(lPtr->Pz) ;
  if ( lPtr->Px != NULL ) free(lPtr->Px) ;
  if ( lPtr->Py != NULL ) free(lPtr->Py) ;

  if ( lPtr->res != NULL ) {
    if ( lPtr->res->monoprog != NULL ) free(lPtr->res->monoprog) ;
    if ( lPtr->res->analprog != NULL ) free(lPtr->res->analprog) ;
    if ( lPtr->res->datasrc != NULL ) free(lPtr->res->datasrc) ;
    free(lPtr->res) ;
  }
}
int Res3_GetConvNpts(Conv3_conv3 *cv)
{
  int i, n ;
  if( ! cv || ! cv->convpts || cv->ndata < 1 ) return 0 ;
  n = 0 ;
  for( i=0 ; i<cv->ndata ; i++ ) n += cv->convpts[i].nres ;
  return n ;
}



static int Copyconvpts(Conv3_convpt *destPtr, Conv3_convpt *srcPtr)
{
  int i ;
  Conv3_respt *respt ;
  if ( srcPtr == NULL || destPtr == NULL ) return 0 ;
  /* save the respt so it doesnt get overwritten */
  respt = destPtr->respts ;
  *destPtr = *srcPtr ;
  destPtr->respts = respt ;
  if( destPtr->nres < 1 ) return 1 ;
  destPtr->respts = 
    (Conv3_respt *)realloc(respt, destPtr->nres*sizeof(Conv3_respt)) ;
  if( ! destPtr->respts ) return 0 ;
  destPtr->nalloc = destPtr->nres ;
  for( i=0 ; i<destPtr->nres ; i++ ) destPtr->respts[i] = srcPtr->respts[i] ;
  return 1 ;
}

int Res3_CopyConvDatapt(Conv3_conv3 *s, int is, Conv3_conv3 *d, int id)
{
  if( ! s || ! d ) return 0 ;
  if( is < 0 || is >= s->ndata ) return 0 ;
  if( id < 0 || id >= d->ndata ) return 0 ;
  d->y[id] = s->y[is] ;
  d->ye[id] = s->ye[is] ;
  d->h[id] = s->h[is] ;
  d->k[id] = s->k[is] ;
  d->l[id] = s->l[is] ;
  d->Ei[id] = s->Ei[is] ;
  d->Ef[id] = s->Ef[is] ;
  d->fs[id] = s->fs[is] ;
  d->omega[id] = s->omega[is] ;
  d->chiA[id] = s->chiA[is] ;
  d->chiB[id] = s->chiB[is] ;
  d->angleA[id] = s->angleA[is] ;
  d->T[id] = s->T[is] ;
  d->H[id] = s->H[is] ;
  d->Pz[id] = s->Pz[is] ;
  d->Px[id] = s->Px[is] ;
  d->Py[id] = s->Py[is] ;
  d->convpts[id].nres = s->convpts[is].nres ;
  if( ! Copyconvpts(d->convpts+id, s->convpts+is) ) return 0 ;
  return 1 ;
}

void Res3_CopyResinfo(Conv3_resinfo *srcRes, Conv3_resinfo *destRes)
{
  if( destRes->monoprog != NULL ) free(destRes->monoprog) ;
  if( destRes->analprog != NULL ) free(destRes->analprog) ;
  if( destRes->datasrc != NULL ) free(destRes->datasrc) ;
  *destRes = *srcRes ;
  if( srcRes->monoprog != NULL ) destRes->monoprog = strdup(srcRes->monoprog) ;
  if( srcRes->analprog != NULL ) destRes->analprog = strdup(srcRes->analprog) ;
  if( srcRes->datasrc != NULL ) destRes->datasrc = strdup(srcRes->datasrc) ;
}
static int
Copyconv3(destPtr, srcPtr)
    Conv3_conv3 *destPtr, *srcPtr;
{
  int i ;
  if ( srcPtr == NULL || destPtr == NULL ) return 0 ;

  if (srcPtr->ndata < 1 || srcPtr->res == NULL) return 0 ;
  if (destPtr->res == NULL) destPtr->res = Res3_NewResinfo() ;
  Res3_CopyResinfo(srcPtr->res, destPtr->res) ;
  if( ! Res3_checkConvNpts(srcPtr->ndata, destPtr) ) return 0 ;
  for( i=0 ; i<srcPtr->ndata ; i++ ) {
    destPtr->y[i] = srcPtr->y[i] ;
    destPtr->ye[i] = srcPtr->ye[i] ;
    destPtr->h[i] = srcPtr->h[i] ;
    destPtr->k[i] = srcPtr->k[i] ;
    destPtr->l[i] = srcPtr->l[i] ;
    destPtr->e[i] = srcPtr->e[i] ;
    destPtr->Ei[i] = srcPtr->Ei[i] ;
    destPtr->Ef[i] = srcPtr->Ef[i] ;
    destPtr->fs[i] = srcPtr->fs[i] ;
    destPtr->omega[i] = srcPtr->omega[i] ;
    destPtr->chiA[i] = srcPtr->chiA[i] ;
    destPtr->chiB[i] = srcPtr->chiB[i] ;
    destPtr->angleA[i] = srcPtr->angleA[i] ;
    destPtr->T[i] = srcPtr->T[i] ;
    destPtr->H[i] = srcPtr->H[i] ;
    destPtr->Px[i] = srcPtr->Px[i] ;
    destPtr->Py[i] = srcPtr->Py[i] ;
    destPtr->Pz[i] = srcPtr->Pz[i] ;
    if( ! Copyconvpts(srcPtr->convpts+i, destPtr->convpts+i) ) return 0 ;
  }
  return 1 ;
}





/*
  computes all eigenvalues and eigenvectors of a real symmetric
  matrix a[0..n-1][0..n-1].  On output elements of a above the diagonal
  are destroyed.  d returns the eigenvalues.
  v is the matrix whose columns contain, 
  on output the normalized eigenvectors of a.
  nrot returns the number of Jacobi rotations required.
  Inverse of matrix a[][] can be obtained as a-1 = U D-1 Uadj
  where U is unitary columns of normalized eigenvectors
  and D-1 is diagonal 1/eigenvalues.
*/

static int jacobi( double **aa, int n, double *d, double **v, int *nrot )
{
  int j, iq, ip, maxitj, nit ;
  int status ;
  double tresh, theta, tau, t, sm, s, h, g, c ;
  double *b, *z ;
  double **a = (double**)0 ;

  a = getmatrix( (double**)0, 0, n, n ) ;

  b = (double *) malloc( (unsigned) n*sizeof(double) ) ;
  z = (double *) malloc( (unsigned) n*sizeof(double) ) ;

  for( ip=0 ; ip<n ; ip++ )
    {
      for( iq=0 ; iq<n ; iq++ ) a[ip][iq] = aa[ip][iq] ;
    }

  for( ip=0 ; ip<n ; ip++ )
    {
      for( iq=0 ; iq<n ; iq++ ) { v[ip][iq] = 0.0 ; }
      v[ip][ip] = 1.0 ;
    }
  for( ip=0 ; ip<n ; ip++ )
    {
      b[ip] = d[ip] = a[ip][ip] ;
      z[ip] = 0.0 ;
    }

  *nrot = 0 ;

  maxitj = 50 ;
  nit = 0 ;
  status = 0 ;

  while( nit < maxitj )
    {
      sm = 0.0 ;
      for( ip=0 ; ip<n-1 ; ip++ )
	{
	  for( iq=ip+1 ; iq<n ; iq++ )
	    {
	      sm += fabs(a[ip][iq]) ;
	    }
	}
      if( sm == 0.0 )  { status = 1 ; break ; }

      if( nit < 3 )
	{
	  tresh = 0.2*sm/(n*n) ;
	}
      else
	{
	  tresh = 0.0 ;
	}


      for( ip=0 ; ip<n-1 ; ip++ )
	{
	  for( iq=ip+1 ; iq<n ; iq++ )
	    {
	      g=100.0*fabs(a[ip][iq]) ;
	      if( nit > 3 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
		          && (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
		{
		  a[ip][iq] = 0.0 ;
		}
	      else if( fabs(a[ip][iq]) > tresh)
		{
		  h = d[iq] - d[ip] ;
		  if( (float)(fabs(h)+g) == (float)fabs(h))
		    {
		      t = (a[ip][iq])/h ;
		    }
		  else
		    {
		      theta = 0.5*h/(a[ip][iq]) ;
		      t = 1.0/(fabs(theta) + sqrt(1.0+theta*theta)) ;
		      if( theta < 0.0 ) t = -t ;
		    }
		  c = 1.0/sqrt(1+t*t) ;
		  s = t*c ;
		  tau = s/(1.0 + c) ;
		  h = t*a[ip][iq] ;
		  z[ip] -= h ;
		  z[iq] += h ;
		  d[ip] -= h ;
		  d[iq] += h ;
		  a[ip][iq] = 0.0 ;

		  for( j=0 ; j<=ip-1 ; j++ )
		    { ROTATE(a,j,ip,j,iq) }
		  for( j=ip+1 ; j<=iq-1 ; j++ )
		    { ROTATE(a,ip,j,j,iq) }
		  for( j=iq+1 ; j<n ; j++ )
		    { ROTATE(a,ip,j,iq,j) }
		  for( j=0 ; j<n ; j++ )
		    { ROTATE(v,j,ip,j,iq) }
		  ++(*nrot) ;
		}
	    }
	}
      ++nit ;
      for( ip=0 ; ip<n ; ip++ )
	{
	  b[ip] += z[ip] ;
	  d[ip]  = b[ip] ;
	  z[ip]  = 0.0   ;
	}
    }
  free(b) ;
  free(z) ;

  getmatrix(a, n, 0, 0) ;
  return status ;
}


/********************************************************************
     RESM calculates the resolution matrix via Cooper and Nathans,
     double-precision version.
     Acta Cryst. 23,357 (1967).

     IFIX < 2 is powder configuration i.e. XMOS(2)=inf.
     If DSP(2).LE.0. or DSP(2).GT.998. OR XMOS(3) is anomalous then
     the analyzer is removed for two-axis configuration. Energy values
     still returned.

     Q is length of scattering vector (A-1).
     EN is neutron energy transfer, positive for neutron energy loss.
     The resolution function is given in terms of the matrix RM(I,J) 
     I=1 Q component along KI-KF 
     I=2 Q component perpendicular to 1 position counterclockwise
     I=3 EN axis and I=4 is vertical Q component.

     HCOL(4) are horizontal collimations FWHM mins from reactor to detector
     VCOL(4) are vertical collimations FWHM degrees from reactor to detector
     MOSA(3) are crystal mosaics FWHM minutes for MONO SAMP and ANAL
     HCO VCO AND MOS are above parameters converted to stnd dev rads.
     Both total norm and norm*volume of ellipsoid are returned
     i.e. multiplied by Jacobian (vol of resolution appropriate for
     integrations in dimensionless variables).
     Normalization as in Chesser and Axe.


  circa 2000
  Modify to handle small crystals, perfect crystals, offset crystals
  This means passing xtal angles, positions and sizes and distances to sample.
  Also return normalized eigenvectors and values, and offet in QE
  xtal-params:  offset-angle-deg   displacements(x  y)
                     length  thickness  distance-to-sample
  just make sure all distances use same units

  xposition	 reactor mono   samp   anal  det    
                     x is along incoming beam to crystal

  yposition      reactor mono   samp   anal  det    
                     y is perp to x so x X y is up.

  height         reactor mono   samp   anal  det    
                     reactor and detector irrelevant

  width          reactor mono   samp   anal  det    
                     reactor-line-source size similar det
                     width along omega, thick perp

  omegaoffset    reactor mono   samp   anal  det    
                     in deg diff from nominal for mono and anal
                     but abs for samp since cant calc nominal here.

  distance       RM  MS  SA  AD	 

  all distances, positions in same units

  NB  most i/o to resm is thru this file globals,
  so this function should be static

  For purposes of doing CS integrals somewhere else,
  save any linear relationships in C&N coordinates
  as they produce delta functions of the same type used to represent the CS's.
  This means that the full uncorrected 3D res ellipsoid 
  should be returned also.
  However, we also want the eigenvectors/values 
  for the corrected resolution function for plotting.
  For this purpose we can't apply the linear constraints as a delta function 
  since we don't want
  to integrate over the Qx Qy E variables.
  Instead we have to find the eigenvectors.

******************************************************************************/


static int resm( double eI, double eF, double q )
{
	static double hco[4], vco[4], mos[3], vmo[3] ;
	static double vcosq[4] ;

	static double ki, kf ;
	static double kisql, kfsql, qsq ;
	static double energ, aom ;

	static double twothet[5], omegal[2], tiltl[2], chit[2] ;
	static double xhco[4], xhcosq[4], xvcosq[4], hcosq[4] ;
	static double hcgsq[4], vcgsq[4] ;
	static double xvco[4], wgeo[4] ;
	static double gammaoffr[4], gammaoffsq[4] ;
	static double deltaoffr[4], deltaoffsq[4], don[4], don2[4] ;
	//static double geoalphar[4], geoalphasq[4] ;
	//static double geobetar[4], geobetasq[4] ;
	static double omM, omM2, omA, omA2, omMv, omAv ;
	static double gon[4], gtn[4], gon2[4] ;
	static double aIn[3], aOut[3], mIn[3], mOut[3] ;
	static double mNorm[3], aNorm[3], omg[4], vmot[2] ;
	static double kifac, kffac ;
	static double deltaoma, deltaomabar, deltaomm, deltaommbar ;
	static double ngg, ng, do2, C, O ;
	static double Bgbar, Bkbar, Aggt, Agkt, Akkt ;
	static double Bd0, Bd1, Bd01, CdX, abet0sq, abet1sq ;
	static double temp ;
	/* static double betam0sqinv, betam1sqinv */
	//static double yrange[5] ;

	static double snm, sna, tom, tomk, toa, toak ;
	static double sinaomMv, sinaomAv ;
	static double lambda, lamsq, beta, alpha ;
	static double bigB, bigA, sb, sa, c ;
	static double e ;
	static double mosq[3], amos[3], amosq[3] ;


	static double c00 ;

	static double vfacM ;
	static double ktsq, rat, fluxperk, fhm ;

	static double a1, a2, a3, a4, a5, a6, a7, a8, a9, aA ;
	static double a3t, a4t, a7t, a8t, a9t, aAt ;
	//static double a11, a12, a22, a33, a44 ;
	//static double a55, a66, a56, a77, a88, a78, a99, aAA, a9A ;
	static double b0, b1, b2, b3, b4, b5 ;
	static double bx1, by1, bx2, by2, bc ;
	static double BM, BA ;
	static double ap ;
	static double bcle[3] ;

	static double s, ss, sv, den, fac ;


	int i, j, k, l ;
	int i1, i2 ;

	int powder = 0 ;
	int anal = 1 ;

	static double pi = 3.14159265 ;
	static double pi2 = 1.570796325 ;
	static double twopi = 6.2831853 ;
	static double rt2pi = 2.5066283  ;	/* sqrt(2.*PI) */
	//static double rtln2 = 0.8325545 ;       /* to convert eigs to HW */
	//static double sqrt2 = 1.4142136 ;
	static double QSMALL = 1.e-6 ;

	static double f  = 4.144283578;	/* hbar**2/massneutron (mev A**2) */
	static double finv = 0.241296229 ;
	static double Dn = 2.072141789 ;
	static double convm = .0001235288 ;	/* FWHM min to sigma radians */
	static double convd = .0074117309 ;	/* FWHM deg to sigma radians */
	//static double convs = .424661 ;
	/* FWHM rad to radians sigma = fwhm/(2*sqrt(2ln2)) */
	static double fwhmtosig = 0.4246609 ;
	static double degtorad = .017453293 ;
	//static double mintorad = .000290888 ;
	/*     CONVM = PI/180/60/(2SQRT(2LN2))  FWHM MIN TO SIGMA RADIANS */
	/*     CONVD = 60 * CONVM               FWHM DEG TO SIGMA RADIANS */


	double vmom, vmomsq, vmoa, vmoasq ;
	double b0msq, b3asq, lmvsq, lavsq ;
	double rmzz ;
	double tiltm, tilta ;
	double dv00, dv30, d0, d3 ;
	double diinv, djinv, ntm1, ntm0, nta3, nta2 ;

	double xcomp, ycomp, zcomp ;
	double dist, s2thet, c2thet ;

	double adiag[3][3] ;

	//static double **adiag = (double**)0 ;

	/* one time allocation of ancilliary matrices */
	//if( adiag == (double**)0 )
	//{
	//  adiag = getmatrix((double**)0, 0, 3, 3) ;
	//}

	/* convert to sigma units from FWHM & convert zero coll to infinity */

	for( i=0 ; i<4 ; ++i )
	  {
	    if( hcol[i] <= 0. ) 	hco[i] = 9999. ;
	    else			hco[i] = convm*hcol[i] ;
	    hcosq[i] = hco[i]*hco[i] ;
	    if( vcol[i] <= 0. )	vco[i] = 9999. ;
	    else			vco[i] = convd*vcol[i] ;
	    vcosq[i] = vco[i]*vco[i] ;
	  }

	for( i=0 ; i<3 ; ++i )
	  {

	    if( i != 1 && mosa[i] <= 0. )
	      mos[i] = 1.e-2 * convm ; /*arc sec def*/
	    else 	    mos[i] = convm*mosa[i] ;
	    mosq[i] = mos[i]*mos[i] ;
	    amos[i] = 1./mos[i] ;
	    amosq[i] = amos[i]*amos[i] ;
	    if( vmos[i] <= 0. ) vmo[i] = mos[i] ;
	    else vmo[i] = convm*vmos[i] ;
	    for( j=0 ; j<3 ; j++ ) rm[i][j] = 0. ;
	  }

	if( eI <= 0. || eF <= 0. )
	  {
	    setErrmsg("ERROR: resm, neutron energies must be > zero!") ;
	    return (0) ;
	  }

	if( q < QSMALL )
	  {
	    setErrmsg("ERROR: resm, Can't calc resolution for Q=0!") ;
	    return (0) ;
	  }
	/*
	  qmin = fabs(ki-kf) ; qmax = ki + kf ;
	  if( q < qmin || q > qmax )
	  {
	  strcat(errmsg," !ERR Can't close scattering triangle!") ;
	  return (0) ;
	  }
	*/

	/* check for powder option */

	powder = 0 ;
	energ = eI - eF ; qsq = q*q ;
	if( mos[1] < 0. || mos[1] > 997. )
		{ energ=0. ; powder = 1 ; }

	/* check for analyzer */

	anal = 1 ;
	if( dsps[1]<=0. || dsps[1]>998. || mosa[2]>998. ) {
	  energ=0. ;
	  anal = 0 ;
	  eF = eI ;
	  mos[2] = 1.e12 ; mosq[2] = mos[2]*mos[2] ;
	  amosq[2] = 0. ; amos[2] = 0. ;
	} else {
	  amos[2] = 1./mos[2] ; amosq[2] = amos[2]*amos[2] ;
	}
	if( dsps[0] <= 0. ) {
	  setErrmsg("ERROR: resm dsp-mono must be greater than zero!") ;
	  return (0) ;
	}

	ki = sqrt(eI/Dn) ;
	kf = sqrt(eF/Dn) ;

	snm = spin[0]*pi/(dsps[0]*ki) ;		 /* sin(THETAmono) */

	if( fabs(snm) > 1. ) {
	  setErrmsg("ERROR: resm, Can't get kI with this dsp-mono!") ;
	  return (0) ;
	}


	twothet[1] = 2.*asin(snm) ;
	tom = snm/sqrt(1.-snm*snm) ;	        /* tan(THETAmono) */


	finv = 1./f ;
	if( anal ) sna = spin[2]*pi/(dsps[1]*kf) ;	/* sin(THETAanal) */
	else sna = 0. ;
	if( fabs(sna) >= 1. ) {
	  if( anal ) {
	    setErrmsg("ERROR: resm, Can't get kF with this dsp-anal!") ;
	    return (0) ;
	  }
	  sna = 0.00001 ;
	}
	twothet[3] = 2.*asin(sna) ;
	toa = sna/sqrt(1.-sna*sna) ;	        /* tan(THETAanal) */

	if( ! anal ) vmo[2] = 0. ;


	aom = 2.*energ*finv ;

	kisql = ki*ki ;  kfsql = kf*kf ;

	/*
	  original C&N is for LRL.  To do L at sample equivalent is to
	  change sign of all sin() terms above.  This still gives Qy
	  as rhc wrt to Qx and Qz up
	  spin 1 -1 1 is supposed to be lrl C&N deflt case
	  all angles positive is ccw from above, so left is pos
	  phi is the acute angle between ki and -q
	*/


	beta = -(qsq-2.*kisql+aom)/(2*ki*kf) ;	/* cos(2thetsam) */
	if( fabs(beta) >= 1. )
	  {
	    setErrmsg("ERROR: resm, invalid scattering triangle!") ;
	    return (0) ;
	  }
	alpha = spin[1]*sqrt(1.-beta*beta) ;	/* sin(2thetsam) */
	twothet[2] = asin(alpha) ;
	bigB  = -(qsq-aom)/(2.*q*kf) ;		/* cos(2thetsam+phi) */
	if( fabs(bigB) >= 1. )
	  {
	    setErrmsg("ERROR: resm, invalid scattering triangle!") ;
	    return (0) ;
	  }
	bigA  = spin[1]*sqrt(1.-bigB*bigB) ;	/* sin(2thetsam+phi) */
	sb = (qsq+aom)/(2.*q*ki) ;		/* cos(phi) */
	if( fabs(sb) >= 1. )
	  {
	    setErrmsg("ERROR: resm, invalid scattering triangle!") ;
	    return (0) ;
	  }
	sa = spin[1]*sqrt(1.-sb*sb) ;		/* sin(phi) */

	/*
	  calc geometrical offsets for deviation angles 
	  and xtal-size effective collimation
	  also combine Soller and geo collimation
	*/

	twothet[0] = pi2 ; twothet[4] = pi2 ;
	/* reactor and detector at zero, use length for yrange */

	/* omegaoff is omegaNominalbl - omegaActualbl = OmegamOff wrt beam ln */
	for( i=0 ; i<2 ; i++ ) {
	  omegal[i] = degtorad*omegaoff[i] ;
	  tiltl[i] = degtorad*tilts[i] ;
	}

	xposition[0] = 0. ; /* source displacement along RM taken zero */
	yposition[0] = 0. ;
	zposition[0] = 0. ;
	xposition[4] = 0. ; /* det displace along AD taken zero */

	/*
	  no longitudinal offset position for reactor and detector
	  distance[0] reactor to mono ...
	  xyposition[0] source in mono coords already
	  r[i](i+1 coords) = riy cos(2theti) - rix sin(2theti)
	  angles positive ccw
	  xihat along path i-1 to i
	  yihat perp ccw
	  angle 1 = mono  angle 2 = sample etc
	*/

	/*
	for( i=0 ; i<5 ; i++ ) {
	  angle = twothet[i]/2. - omegal[i] ;
	  yrange[i] = length[i]*fabs(sin(angle)) + width[i]*cos(angle) ;
	}
	*/

	/*
	  April 2003
	  Change global inputs to gammaoff deltaoff geoalpha geobeta
	  in FWHM degrees

	  Dec 2009
	  there are 2 basic offset-input methods for crystal offsets:
	  angle-mode and position mode.
	  the default is angle mode.
	  if no flags are set for offset method, the default is to
	  input all the angles in degrees
	  gammaoffsets and deltaoffsets for crystals
	  omegaoff = omegaNom - omegaAct wrt beam line
	  tilts
	  Soller and geometrical collimations.
	  path angle-offsets after a crystal are related to the
	  path angle-offsets before the crystal and the crystal angle offsets
	  This constraint is always enforced.
	  In angle-mode there is no way to check that the beam is aimed
	  at the sample, which or course it must for the calculations
	  made here. 
	  On the other hand the anlyzed beam does not
	  have to aim at a detector point. This could handle a PSD as
	  long as pixel efficiency corrections are made to the data.

	  Position mode uses offset positions and
	  distances between elements to calculate the angles.
	  reactor offset position is constrained to zero, but
	  detector offset can be nonzero to handle PSD case.

	  If lengths and widths are non-zero the corresponding
	  geometrical collimation is calculated if geomode flag is set
	*/


	for( i=0 ; i<4 ; i++ )
	  {
	    /* convert degrees to rad */
	    gammaoffr[i] = degtorad*gammaoff[i] ;
	    deltaoffr[i] = degtorad*deltaoff[i] ;
	    /* gammaoffsq[i] = gammaoffr[i]*gammaoffr[i] */
	    /* deltaoffsq[i] = deltaoffr[i]*deltaoffr[i] */
	  }



	/*
	  gammaoff10 is related to gammaoff00
	  and same for gammaoff30 20
	  same for deltaoffs
	  etc
	  enforce this condition in angle input mode
	*/

	/* effective xtal setting angles */
	omM = twothet[1]/2. ;
	omA = twothet[3]/2. ;

	if( offseton ) {
	  omM = twothet[1]/2. - omegal[0] - gammaoffr[0] ;
	  omA = twothet[3]/2. - omegal[1] - gammaoffr[2] ;
	}


	/* make the vertical corrections to omega effective */
	sinaomMv = sin(fabs(omM)) ;
	if( offseton ) sinaomMv = cos(deltaoffr[0])*cos(tiltl[0])*sin(fabs(omM))
	  -sin(deltaoffr[0])*sin(tiltl[0]) ;
	if( omM < 0 ) omMv = -asin(sinaomMv) ;
	else omMv = asin(sinaomMv) ;
	kie = pi/(dsps[0]*sinaomMv) ;
	tom = tan(omMv) ;

	if( anal ) {
	  sinaomAv = sin(fabs(omA)) ;
	  if( offseton ) sinaomAv =
			   cos(deltaoffr[2])*cos(tiltl[1])*sin(fabs(omA))
			   -sin(deltaoffr[2])*sin(tiltl[1]) ;
	  if( omA < 0 ) omAv = -asin(sinaomAv) ;
	  else omAv = asin(sinaomAv) ;
	  kfe = pi/(dsps[1]*sinaomAv) ;
	  toa = tan(omAv) ;
	} else {
	  sinaomAv = 0. ;
	  kfe = kf ;
	  toa = 0. ;
	}

	/* enforce relations between pre and post offset path angles */
	gammaoffr[1] = -gammaoffr[0] - 2.*omegal[0] ;
	gammaoffr[3] = -gammaoffr[2] - 2.*omegal[1] ;
	/* use small angle approx */
	deltaoffr[1] = deltaoffr[0] + 2.*sin(fabs(omMv))*sin(tiltl[0]) ;
	deltaoffr[3] = deltaoffr[2] + 2.*sin(fabs(omAv))*sin(tiltl[1]) ;



	/*
	  monomode analmode
	  mode 0   default uses input angles
	  mode 1   use positions to calc offsets
	*/

	if( offseton && monomode && distance[0] > 0. && distance[1] > 0. ) {
	  /* use mono xtal position to calculate angles */
	  /* compute nominal vector onto mono in mono coords */
	  /* normally xposition[0]=yposition[0]=0 ie source not displaced */
	  mIn[0] = distance[0] + xposition[1] - xposition[0] ;
	  mIn[1] = yposition[1] - yposition[0] ;
	  mIn[2] = zposition[1] - zposition[0] ;
	  dist = unitvec(mIn) ;
	  deltaoffr[0] = asin(mIn[2]) ;
	  gammaoffr[0] = asin(mIn[1]/cos(deltaoffr[0])) ;

	  /* compute nominal vector onto sample in mono coords */
	  /* but need sample coord to get post xtal gammaoff */
	  /* [xyz]position[2] = 0 ie sample not displaced */
	  s2thet = sin(twothet[1]) ;
	  c2thet = cos(twothet[1]) ;
	  mOut[0] = distance[1]*c2thet - xposition[1] ;
	  mOut[1] = distance[1]*s2thet - yposition[1] ;
	  mOut[2] = -zposition[1] ;
	  dist = unitvec(mOut) ;
	  if( dist <= 0. ) dist = 1. ;

	  deltaoffr[1] = asin(mOut[2]) ;
	  /* can get gammaoffr[1] from yhatS(inMonoCoord) . mOut */
	  /* yhatS = cos2tm yhatM - sin2tm xhatM */
	  ycomp = s2thet*xposition[1] - c2thet*yposition[1] ;
	  gammaoffr[1] = asin(ycomp/(dist*cos(deltaoffr[1]))) ;
	  
	  /* get the xtal norm unit vec * 2sin(omeff) from mOut - mIn */
	  for( i=0 ; i<3 ; i++ ) mNorm[i] = mOut[i] - mIn[i] ;
	  dist = unitvec(mNorm) ;
	  omMv = asin(dist/2.) ;
	  sinaomMv = sin(omMv) ;
	  if( twothet[1] < 0. ) omMv = -omMv ;
	  kie = pi/(dsps[0]*sinaomMv) ;
	  tom = tan(omMv) ;
	  tiltl[0] = mNorm[2] ;
	  omegal[0] = -(gammaoffr[0] + gammaoffr[1])/2. ;
	  omM = twothet[1]/2. - omegal[0] - gammaoffr[0] ;
	}

	if( offseton && anal
	    && analmode && distance[2] > 0. && distance[3] > 0. ) {
	  /* use analyzer xtal position to calc angles */
	  /* compute nominal vector onto analyzer in analyzer coords */
	  /* normally [xyz]position[2] = 0 ie sample not displaced */
	  aIn[0] = distance[2] + xposition[3] ;
	  aIn[1] = yposition[3] ;
	  aIn[2] = zposition[3] ;
	  dist = unitvec(aIn) ;
	  deltaoffr[2] = asin(aIn[2]) ;
	  gammaoffr[2] = asin(aIn[1]/cos(deltaoffr[2])) ;

	  /* compute nominal vector onto detector in analyzer coords */
	  /* but need det coord to get post xtal gammaoff */
	  /* normally xposition[4]=0 but [yz]position[4] can be non-zero PSD */
	  s2thet = sin(twothet[3]) ;
	  c2thet = cos(twothet[3]) ;
	  aOut[0] = distance[3]*c2thet - yposition[4]*s2thet - xposition[3] ;
	  aOut[1] = distance[3]*s2thet + yposition[4]*c2thet - yposition[3] ;
	  aOut[2] = zposition[4] - zposition[3] ;
	  dist = unitvec(aOut) ;
	  if( dist <= 0. ) dist = 1. ;

	  deltaoffr[3] = asin(aOut[2]) ;
	  /* get gammaoffr[3] from yhatD . aOut */
	  /* yHatD = cos2tA yhatA - sin2tA xhatA */
	  /* ycomp = -s2thet*
	    (distance[3]*c2thet - yposition[4]*s2thet - xposition[3])
	    + c2thet*
	    (distance[3]*s2thet + yposition[4]*c2thet - yposition[3]) ; */
	  ycomp = yposition[4] + s2thet*xposition[3] - c2thet*yposition[3] ;
	  gammaoffr[1] = asin(ycomp/(dist*cos(deltaoffr[3]))) ;
	  
	  /* get the xtal norm unit vec * 2sin(omeff) from aOut - aIn */
	  for( i=0 ; i<3 ; i++ ) aNorm[i] = aOut[i] - aIn[i] ;
	  dist = unitvec(aNorm) ;
	  omAv = asin(dist/2.) ;
	  sinaomAv = sin(omAv) ;
	  kfe = pi/(dsps[1]*sinaomAv) ;
	  toa = tan(omAv) ;
	  if( twothet[3] < 0. ) omAv = -omAv ;
	  tiltl[1] = aNorm[2] ;
	  omegal[1] = -(gammaoffr[2] + gammaoffr[3])/2. ;
	  omA = twothet[3]/2. - omegal[1] - gammaoffr[2] ;
	}


	/* compute the effective crystal setting angles omegaEffm and a */
	/* these are needed for calc of geometrical collimation */


	for( i=0 ; i<4 ; i++ )
	  {
	    /* convert from FWHM degrees to sigma rad */
	    xhco[i] = convd*geoalpha[i] ;
	    xvco[i] = convd*geobeta[i] ;
	    gammaoffsq[i] = gammaoffr[i]*gammaoffr[i] ;
	    deltaoffsq[i] = deltaoffr[i]*deltaoffr[i] ;
	  }


	if( geomode ) {
	  /* calc the geo collimation */
	  omg[0] = fabs(omM) ;
	  omg[1] = fabs(twothet[1]/2. + omegal[0] + gammaoffr[1]) ;
	  omg[2] = fabs(omA) ;
	  omg[3] = fabs(twothet[3]/2. + omegal[1] + gammaoffr[3]) ;
	  
	  for( i=0 ; i<4 ; i++ ) {
	    if( distance[i] <= 0. ) continue ;
	    xhco[i] = fwhmtosig*(width[2*((i+1)/2)] +
			   width[1+2*(i/2)]*sin(omg[i]) +
			   thick[1+2*(i/2)]*cos(omg[i]))/(2*distance[i]) ;
	    xvco[i] = fwhmtosig*(height[i] + height[i+1])/(2*distance[i]) ;
	  }
	}

	for( i=0 ; i<4 ; i++ )
	  {
	    if( xhco[i] <= 0. ) xhco[i] = 1.e12 ; /* removes it */
	    if( xvco[i] <= 0. ) xvco[i] = 1.e12 ; /* removes it */
	    xhcosq[i] = xhco[i]*xhco[i] ;
	    xvcosq[i] = xvco[i]*xvco[i] ;

	    /* compute Soller-geometrical averages */
	    hcgsq[i] = hcosq[i]/(1. + hcosq[i]/xhcosq[i]) ;
	    vcgsq[i] = vcosq[i]/(1. + vcosq[i]/xvcosq[i]) ;
	  }


	/*
	  compute the vertical offset correction to the crystal setting angles
	  Note that these are second order corrections
	*/



	chit[0] = sinaomMv*tiltl[0] ;
	chit[1] = sinaomAv*tiltl[1] ;

	vmot[0] = sinaomMv*vmo[0] ;
	vmot[1] = sinaomAv*vmo[2] ;


	if( ! anal ) {
	  omegal[1] = 0. ;
	  tiltl[1] = 0. ;
	  chit[1] = 0. ;
	  for( i=2 ; i<=3 ; i++ ) {
	    gammaoffr[i] = 0. ;
	    gammaoffsq[i] = 0. ;
	    deltaoffr[i] = 0. ;
	    deltaoffsq[i] = 0. ;
	  }
	}

	lambda = ki/kf ;			/* lambda C&N */
	lamsq = lambda*lambda ;
	c  = (lambda-beta)/alpha ;    	/* -C C&N as alpha is signed  my mu1 */
	e  = (beta*lambda-1.)/alpha ;   /* -E C&N as alpha is signed  my mu2 */

	/* at this point we can correct ki and kf to kie kfe */
	kinom = ki ;
	kfnom = kf ;
	Ki = kinom ;
	Kf = kfnom ;
	ki = kie ;
	kf = kfe ;
	kisql = ki*ki ;  kfsql = kf*kf ;
	kioff = kie - kinom ;
	kfoff = kfe - kfnom ;

	/* first order corrections to Q,E C&N coords */
	/*
	  qzcnoff = kfe*deltaoffr[2] - kie*deltaoffr[1] ;
	  xyecnoff[0] = -bigA*kfe*gammaoffr[2] + sa*kie*gammaoffr[1]
	  + bigB*kfoff - sb*kioff ;
	  xyecnoff[1] = bigB*kfe*gammaoffr[2] - sb*kie*gammaoffr[1]
	  + bigA*kfoff - sa*kioff ;
	  xyecnoff[2] = Dn*(kioff*(kie+kinom) - kfoff*(kfe+kfnom)) ;
	*/

	toak = toa/kf ;
	tomk = tom/ki ;
	/*
	  toaksq = toak*toak ;
	  tomksq = tomk*tomk ;
	*/

	/*
	  compute aij coef products required as in C&N
	  note signs reverse between a5 and a6 and also a9 and a10
	  this doesnt change C&N result
	  and makes it easier to put in omega offsets
	*/

	a1 = tomk*amos[0] ;
	a2 = amos[0]/ki ;
	a3 = 1./ki/hco[1] ;
	a3t= 1./ki/xhco[1] ;
	a4 = 1./kf/hco[2] ;
	a4t= 1./kf/xhco[2] ;
	a5 = -toak*amos[2] ;
	a6 = amos[2]/kf ;
	a7 = 2.*tomk/hco[0] ;
	a7t= 2.*tomk/xhco[0] ;
	a8 = 1./ki/hco[0] ;
	a8t= 1./ki/xhco[0] ;
	a9 = -2.*toak/hco[3] ;
	a9t= -2.*toak/xhco[3] ;
	aA = 1./kf/hco[3] ;
	aAt= 1./kf/xhco[3] ;


	/*
	  (kI/kIe - 1)(omm - omM) is typically second order effect
	  deltaomegam - gamma10 = omegal[0]
	*/
	kifac = (Ki/kie - 1.) ;
	deltaomm = omM - twothet[1]/2. ;
	deltaommbar = kifac*deltaomm ;
	gtn[1] = (omegal[0] + deltaommbar)*amos[0] ; /* omMhat */
	gtn[0] = 2*deltaommbar/hco[0] ; /* om0hat */

	kffac = (Kf/kfe - 1.) ;
	deltaoma = omA - twothet[3]/2. ;
	deltaomabar = kffac*deltaoma ;
	gtn[2] = (omegal[1] - deltaomabar)*amos[2] ; /* omAhat */
	gtn[3] = -2*deltaomabar/hco[3] ; /* om3hat */

	omM2 = gtn[1]*amos[0] ;
	omA2 = gtn[2]*amos[2] ;
	for( i=0 ; i<4 ; i++ ) {
	  gon[i] = gammaoffr[i]/xhco[i] ;
	  gon2[i] = gon[i]/xhco[i] ;
	  don[i] = deltaoffr[i]/xvco[i] ;
	  don2[i] = don[i]/xvco[i] ;
	}
	gon[0] -= 2*deltaommbar/xhco[0] ;
	gon[3] += 2*deltaomabar/xhco[3] ;

	/*
	  a5,6,9,A change sign wrt C&N but they occur in the b's
	  below as aAA a66 a55 a99 a56 a9A
	  so the b's are unaffected

	a11 = tomksq*amosq[0] ;
	a22 = amosq[0]/kisql ;
	a12 = tomk*amosq[0]/ki ;
	a33 = 1./kisql/hcgsq[1] ;
	a44 = 1./kfsql/hcgsq[2] ;
	a55 = toaksq*amosq[2] ;
	a66 = amosq[2]/kfsql ;
	a56 = -toak*amosq[2]/kf ;
	a77 = 4.*tomksq/hcgsq[0] ;
	a88 = 1./kisql/hcgsq[0] ;
	a78 = 2.*tomk/ki/hcgsq[0] ;
	a99 = 4.*toaksq/hcgsq[3] ;
	aAA = 1./kfsql/hcgsq[3] ;
	a9A = -2.*toak/kf/hcgsq[3] ;
	*/
	
	b0 = a1*a2 + a7*a8 + a7t*a8t ;
	b5 = a1*a1 + a7*a7 + a7t*a7t ;

	b1 = a2*a2 + a3*a3 + a8*a8 + a3t*a3t + a8t*a8t ;
	b2 = a4*a4 + a6*a6 + aA*aA + a4t*a4t + aAt*aAt ;
	b3 = a5*a5 + a9*a9 + a9t*a9t ; /* zero no anal */
	b4 = a5*a6 + a9*aA + a9t*aAt ; /* zero no anal */

	c2bar[0][1] = 0. ;
	c2bar[0][2] = 0. ;
	c2bar[1][0] = 0. ;
	c2bar[2][0] = 0. ;
	c2bar[0][0] = b1 ;
	c2bar[1][1] = b2 ;
	c2bar[2][2] = b3 ;
	c2bar[1][2] = b4 ;
	c2bar[2][1] = b4 ;

	/* we also need new linear terms for Q offset */
	
	bx1 = a1*gtn[1] + a7*gtn[0] - a7t*gon[0] ; /* Kx */
	by1 = a2*gtn[1] + a8*gtn[0] - a8t*gon[0] - a3t*gon[1] ; /* Ky */
	bx2 = a5*gtn[2] + a9*gtn[3] - a9t*gon[3] ; /* Lx */
	by2 = a6*gtn[2] + aA*gtn[3] - aAt*gon[3] - a4t*gon[2] ; /* Ly */
	bcle[0] = by1 ;
	bcle[1] = by2 ;
	bcle[2] = bx2 ;
	bc = bx1 + by1*c + bx2*lambda + by2*e ; /* B0bar */

	/* compute C&N A' */
	ap = b5 + 2*b0*c + b1*c*c + b3*lamsq + 2.*b4*lambda*e + b2*e*e ;

	/* compute the terms used in the g coefs */
	bbar[0] = b0 + b1*c ;
	bbar[1] = b2*e + b4*lambda ;
	bbar[2] = b3*lambda + b4*e ;
	//b0c = b0 + b1*c ;
	//bel = b2*e + b4*lambda ;
	//ble = b3*lambda + b4*e ; /* zero for noanal */

	for( i=0 ; i<3 ; i++ ) {
	  /* also the linear D F H coefs are need to calc the dQE offset */
	  gl[i] = bcle[i] - bc*bbar[i]/ap ;
	  /* gl[2] = 0 no anal */
	  for( j=0 ; j<3 ; j++ ) {
	    gc[i][j] = c2bar[i][j] - bbar[i]*bbar[j]/ap ;
	  }
	}
	/*
	  gc00 = C&N g0
	  gc11 = C&N g1
	  gc22 = C&N g2   zero no anal
	  gc01 = C&N g4/2
	  gc02 = C&N g5/2 zero no anal
	  gc12 = C&N g3/2 zero no anal

	  gc[0][0] = b1 - bbar[0]*bbar[0]/ap ;
	  gc[1][1] = b2 - bbar[1]*bbar[1]/ap ;
	  gc[2][2] = b3 - bbar[2]*bbar[2]/ap ;
	  gc[0][1] = -bbar[0]*bbar[1]/ap ;
	  gc[1][0] = gc[0][1] ;
	  gc[0][2] = -bbar[0]*bbar[2]/ap ;
	  gc[2][0] = gc[0][2] ;
	  gc[1][2] = b4 - bbar[1]*bbar[2]/ap ;
	  gc[2][1] = gc[1][2] ;
	*/


	/*
	  setup the transform matrix DFH = T * dQE
	*/

	dfh[0][0] = -bigB/alpha ;          /* C&N  -d1 */
	dfh[0][1] = -bigA/alpha ;         /* C&N  d2 */
	dfh[0][2] = -finv/alpha/kfnom ;       /* C&N  -d4/hbar */
	dfh[1][0] = -sb/alpha ;            /* C&N  -f1 */
	dfh[1][1] = -sa/alpha ;           /* C&N  f2 */
	dfh[1][2] = -finv*beta/alpha/kfnom ;  /* C&N  -f4/hbar */
	dfh[2][0] = 0. ;
	dfh[2][1] = 0. ;
	dfh[2][2] = -finv/kfnom ;            /* C&N  h4/hbar */

	/*
	  compute the rm as dfhAdj gc dfh
	  NB this is still for the form exp(-1/2 x rm x)
	*/

	for( i=0 ; i<3 ; i++ ) {
	  for( j=0 ; j<=i ; j++ ) {
	    rm[i][j] = 0. ;
	    for( k=0 ; k<3 ; k++ ) {
	      for( l=0 ; l<3 ; l++ ) {
		rm[i][j] += gc[k][l]*dfh[k][i]*dfh[l][j] ;
	      }
	    }
	    rm[j][i] = rm[i][j] ;
	  }
	}

	/* also compute the linear part */
	qzoff = 0. ;
	for( i=0 ; i<3 ; i++ ) {
	  brm[i] = 0 ;
	  xyeoff[i] = 0. ;
	  for( j=0 ; j<3 ; j++ ) {
	    brm[i] += gl[j]*dfh[j][i] ;
	    rminv[i][j] = 0. ;
	  }
	}

	/*
	  the dQE offset vector is -rminv b/2
	  that is the rm origin is at Q0,E0 + dQE
	  so we need rm inverse
	  for 3x3 or 2x2 this seems easy analytically
	  and we need the determinant to do this
	  both 2 and 3 axis have
	  det = (b1*b5 - b0*b0)/ap/sin(2thetS)
	  2-axis then multiplies b2 and
	  3-axis multiplies (b2*b3 - b4*b4)*h4^2
	  NB this is det(rm) from exp(-1/2 x rm x)
	*/

	detrm = (b1*b5 - b0*b0)/ap/fabs(alpha)/fabs(alpha) ;
	if( anal ) {
	  detrm *= (b2*b3 - b4*b4)*finv*finv/kfsql ;
	  if( detrm > 0. ) {
	    rminv[0][0] = (rm[1][1]*rm[2][2] - rm[1][2]*rm[2][1])/detrm ;
	    rminv[1][1] = (rm[0][0]*rm[2][2] - rm[0][2]*rm[2][0])/detrm ;
	    rminv[2][2] = (rm[0][0]*rm[1][1] - rm[0][1]*rm[1][0])/detrm ;
	    rminv[0][1] = (rm[0][2]*rm[2][1] - rm[0][1]*rm[2][2])/detrm ;
	    rminv[1][0] = rminv[0][1] ;
	    rminv[0][2] = (rm[0][1]*rm[1][2] - rm[0][2]*rm[1][1])/detrm ;
	    rminv[2][0] = rminv[0][2] ;
	    rminv[1][2] = (rm[0][2]*rm[1][0] - rm[0][0]*rm[1][2])/detrm ;
	    rminv[2][1] = rminv[1][2] ;
	    if( offseton ) {
	      xyeoff[0] =
		-(rminv[0][0]*brm[0]+rminv[0][1]*brm[1]+rminv[0][2]*brm[2]) ;
	      xyeoff[1] =
		-(rminv[1][0]*brm[0]+rminv[1][1]*brm[1]+rminv[1][2]*brm[2]) ;
	      xyeoff[2] =
		-(rminv[2][0]*brm[0]+rminv[2][1]*brm[1]+rminv[2][2]*brm[2]) ;
	    }
	  }
	} else {
	  detrm *= b2 ;
	  if( detrm > 0. ) {
	    rminv[0][0] = rm[1][1]/detrm ;
	    rminv[1][1] = rm[0][0]/detrm ;
	    rminv[0][1] = -rm[0][1]/detrm ;
	    rminv[1][0] = rminv[0][1] ;
	    if( offseton ) {
	      xyeoff[0] = -(rminv[0][0]*brm[0]+rminv[0][1]*brm[1]) ;
	      xyeoff[1] = -(rminv[1][0]*brm[0]+rminv[1][1]*brm[1]) ;
	      xyeoff[2] = 0. ;
	    }
	  }
	}
	

	/* compute the vertical mosaic and tilt factors */
	vmom = 2.*vmot[0] ;
	vmomsq = vmom*vmom ;
	vmoa = 2.*vmot[1] ; /* zero if no anal */
	vmoasq = vmoa*vmoa ;
	b0msq = vmomsq + vcgsq[0] ;
	b3asq = vmoasq + vcgsq[3] ;
	lmvsq = (1./b0msq + 1./vcgsq[1])/kisql ; /* AMdbar */
	lavsq = (1./b3asq + 1./vcgsq[2])/kfsql ; /* AAdbar */
	rmzz = lmvsq*lavsq/(lmvsq + lavsq) ;

	tiltm = 2.*chit[0] ;
	tilta = 2.*chit[1] ;
	dv00 = (vcgsq[0]/xvcosq[0])*deltaoffr[0] + tiltm ;
	dv30 = (vcgsq[3]/xvcosq[3])*deltaoffr[3] - tilta ;
	d0 = dv00/b0msq + deltaoffr[1]/xvcosq[1] ; /* BMd */
	d3 = dv30/b3asq + deltaoffr[2]/xvcosq[2] ; /* BAd */
	/* need to divide by k to get the BMdbar BAdbar */
	qzoff = 0. ;
	if( offseton ) qzoff = d3/lavsq/kf - d0/lmvsq/ki ;

	/* now we can init the norm r0 from the vertical col facts FVI FVF */
	/* r0 is for CS-> S(Q,E) so cancels 2 Dn Kf from Jacobian */
	BM = sqrt(vcgsq[0]*vcgsq[1]/(vcgsq[0] + vcgsq[1] + vmomsq)) ;
	r0 = BM ;
	if( anal ) {
	  BA = sqrt(vcgsq[2]*vcgsq[3]/(vcgsq[2] + vcgsq[3] + vmoasq)) ;
	} else {
	  BA = sqrt(vcgsq[2]*vcgsq[3]/(vcgsq[2] + vcgsq[3])) ;
	}
	r0 *= BA ;
	r0 *= twopi*mref*sqrt(rmzz/ap)/(ki*kf*fabs(alpha)) ;
	/* include the kf/ki factor so that CS->scatteringLaw */
	r0 *= kf/ki ;
	if( anal ) r0 *= aref ;

	/*
	  start collecting the corrections to the normalization
	  with the horizontal offset corrections

	  Dec 2009
	  multixtal normalization reworked. See multixtal.lyx doc

	*/

	if( offseton ) {
	  ngg = amosq[0] + 1./xhcosq[0] + 1./xhcosq[1] ;
	  ng = amosq[0] + 2./xhcosq[0] ;
	  do2 = 2.*deltaomm ;
	  C = do2*do2/hcosq[0] + gammaoffr[1]*gammaoffr[1]*ngg ;
	  Bgbar = do2/hcosq[0] + gammaoffr[1]*ngg ;
	  Bkbar = 2.*do2/hcosq[0] + gammaoffr[1]*ng ;
	  Aggt = hcgsq[0]*hcgsq[1] ;
	  temp = mosq[0]*hcgsq[1] ;
	  Akkt = Aggt + 4.*temp ;
	  Agkt = Aggt + 2.*temp ;
	  Aggt += temp + mosq[0]*hcgsq[0] ;
	  temp = hcgsq[0] + hcgsq[1] + 4*mosq[0] ;
	  O = Aggt*Bkbar*Bkbar - 2.*Agkt*Bgbar*Bkbar + Akkt*Bgbar*Bgbar;
	  O = (O/temp - C)/2. ;
	  c00 = exp(O) ;
	  
	  if( anal ) {
	    ngg = amosq[2] + 1./xhcosq[2] + 1./xhcosq[3] ;
	    ng = amosq[2] + 2./xhcosq[3] ;
	    do2 = 2.*deltaoma ;
	    C = do2*do2/hcosq[3] + gammaoffr[2]*gammaoffr[2]*ngg ;
	    Bgbar = -do2/hcosq[3] + gammaoffr[2]*ngg ;
	    Bkbar = -2.*do2/hcosq[0] + gammaoffr[2]*ng ;
	    Aggt = hcgsq[2]*hcgsq[3] ;
	    temp = mosq[2]*hcgsq[2] ;
	    Akkt = Aggt + 4.*temp ;
	    Agkt = Aggt + 2.*temp ;
	    Aggt += temp + mosq[2]*hcgsq[3] ;
	    temp = hcgsq[2] + hcgsq[3] + 4*mosq[2] ;
	    O = Aggt*Bkbar*Bkbar - 2.*Agkt*Bgbar*Bkbar + Akkt*Bgbar*Bgbar;
	    O = (O/temp - C)/2. ;
	    c00 *= exp(O) ;
	  }
	  
	  /* now vertical offset correction factors */
	  Bd0 = 4*chit[0]/vmomsq ;
	  CdX = Bd0*chit[0] + deltaoffsq[0]/xhcosq[0]
	    + deltaoffsq[1]/xhcosq[1] ;
	  Bd1 = -Bd0 -2.*deltaoffr[1]/xhcosq[1] ;
	  Bd0 -= 2.*deltaoffr[0]/xhcosq[0] ;
	  abet0sq = 1./vcgsq[0] + 1./vmomsq ;
	  abet1sq = 1./vcgsq[1] + 1./vmomsq ;
	  Bd01 = Bd0*Bd0*abet0sq + Bd1*Bd1*abet1sq + 2.*Bd0*Bd1/vmomsq ;
	  O = (BM*BM*vmomsq*Bd01/4. - CdX)/2. ;
	  c00 *= exp(O) ;
	  
	  if( anal ) {
	    Bd0 = 4.*chit[1]/vmoasq ;
	    CdX = Bd0*chit[1] + deltaoffsq[3]/xhcosq[3]
	      + deltaoffsq[2]/xhcosq[2];
	    Bd1 = Bd0 -2.*deltaoffr[2]/xhcosq[2] ;
	    Bd0 = -Bd0 - 2.*deltaoffr[3]/xhcosq[3] ;
	    abet0sq = 1./vcgsq[3] + 1./vmoasq ;
	    abet1sq = 1./vcgsq[2] + 1./vmoasq ;
	    Bd01 = Bd0*Bd0*abet0sq + Bd1*Bd1*abet1sq + 2.*Bd0*Bd1/vmoasq ;
	    O = (BA*BA*vmoasq*Bd01/4. - CdX)/2. ;
	    c00 *= exp(O) ;
	  }

	  r0 *= c00 ;
	}

	/*
	for( i=0 ; i<4 ; i++ ) c00 += gon[i]*gon[i] ;
	c00 += gtn[0]*gtn[0] + gtn[1]*gtn[1] ;
	if( anal ) c00 += gtn[2]*gtn[2] + gtn[3]*gtn[3] ;
	*/

	/*
	bmg = gon2[0] + gon2[1] - omM2 ;
	bmk = 2*gon2[0] - omM2 ;
	*/

	/*
	ntm0 = gon2[0] - omM2 ;
	ntm1 = ntm0 + gon2[1] ;
	ntm0 += gon2[0] ;

	diinv = mosq[0]*hcgsq[0]*hcgsq[1]/(4.*mosq[0] + hcgsq[0] + hcgsq[1]) ;
	*/

	/*
	gx = gammaoffr[0]/xhcosq[0] ;
	ntm1 = amosq[0]*omegal[0] - gx ;
	ntm0 = ntm1 - gx ;
	ntm1 -= gammaoffr[1]/xhcosq[1] ;
	fac = (4./hcgsq[0] + amosq[0])*ntm1*ntm1 ;
	fac += (1./hcgsq[0] + 1./hcgsq[1] + amosq[0])*ntm0*ntm0 ;
	fac -= 2.*(2./hcgsq[0] + amosq[0])*ntm0*ntm1 ;
	*/

	/*
	fac = (4./hcgsq[0] + amosq[0])*ntm1*ntm1 ;
	fac += (1./hcgsq[0] + 1./hcgsq[1] + amosq[0])*ntm0*ntm0 ;
	fac -= 2.*(2./hcgsq[0] + amosq[0])*ntm1*ntm0 ;
	c00 -= diinv*fac ;

	if( anal ) {
	  djinv = mosq[2]*hcgsq[2]*hcgsq[3]/(4.*mosq[2] + hcgsq[2] + hcgsq[3]);

	    gx = gammaoffr[3]/xhcosq[3] ;
	    nta3 = -amosq[2]*omegal[1] + gx ;
	    nta2 = -nta3 + gx ;
	    nta3 -= gammaoffr[2]/xhcosq[2] ;

	  nta3 = gon2[3] - omA2 ;
	  nta2 = nta3 + gon2[2] ;
	  nta3 += gon2[3] ;

	  fac = (4./hcgsq[3] + amosq[2])*nta2*nta2 ;
	  fac += (1./hcgsq[3] + 1./hcgsq[2] + amosq[2])*nta3*nta3 ;
	  fac -= 2.*(2./hcgsq[3] + amosq[2])*nta3*nta2 ;
	  c00 -= djinv*fac ;
	}
	*/

	/*
	  bmd0 = omM2 - 2.*gon2[0] ;
	  bmd1 = -omM2 - 2.*gon2[1] ;
	*/

	/*
	  betam0sqinv = 1./vcgsq[0] + 1./vmomsq ;
	  betam1sqinv = 1./vcgsq[1] + 1./vmomsq ;
	*/

	/* now the vertical offset corrections */

	/*
	c00 += 4.*chit[0]*chit[0]/vmomsq + deltaoffr[0]*don2[0]
	  + deltaoffr[1]*don2[1] ;
	diinv = 0.25*vmomsq*vcgsq[0]*vcgsq[1]/(vcgsq[0] + vcgsq[1] + vmomsq) ;
	ntm0 = 4.*chit[0]/vmomsq ;
	ntm1 = -ntm0 - 2.*don2[1] ;
	ntm0 -= 2.*don2[0] ;

	fac = (1./vcgsq[0] + 1./vmomsq)*ntm0*ntm0 ;
	fac += (1./vcgsq[1] + 1./vmomsq)*ntm1*ntm1 ;
	fac += 2.*ntm0*ntm1/vmomsq ;
	c00 -= diinv*fac ;

	if( anal ) {
	  c00 += 4.*chit[1]*chit[1]/vmoasq + deltaoffr[3]*don2[3]
	    + deltaoffr[2]*don2[2] ;
	  djinv = 0.25*vmoasq*vcgsq[3]*vcgsq[3]/(vcgsq[3] + vcgsq[2] + vmoasq);
	  nta2 = 4.*chit[1]/vmoasq ;
	  nta3 = -nta2 - 2.*don2[3] ;
	  nta2 -= 2.*don2[2] ;
	  
	  fac = (1./vcgsq[3] + 1./vmoasq)*nta3*nta3 ;
	  fac += (1./vcgsq[2] + 1./vmoasq)*nta2*nta2 ;
	  fac += 2.*nta3*nta2/vmoasq ;
	  c00 -= djinv*fac ;
	}
	*/

	/*
	  c00 -= d0*d0/lmvsq/kisql + d3*d3/lavsq/kfsql ;
	  c00 += betaoffsq[0]/(vcosq[0] + xvcosq[0]) ;
	  c00 += dv00*dv00/b0msq + betaoffsq[1]/xvcosq[1] ;
	  c00 += dv30*dv30/b3asq + betaoffsq[2]/xvcosq[2] ;
	  c00 += betaoffsq[3]/(vcosq[3] + xvcosq[3]) ;
	*/

	/*
	  We cant treat perfect xtals since intensity goes to zero
	*/

	/* ********************************
	  Now do sample mosaic corrections
	  use the Q value corrected for offset
	*/

	qsq = q+xyeoff[0] ;
	qsq = sqrt(qsq*qsq + xyeoff[1]*xyeoff[1]) ;
	s = qsq*mos[1]*mos[1] ;
	sv = qsq*vmo[1]*vmo[1] ;
	den = s*rm[1][1] ;

	if( powder )
	  {
	    if( rm[1][1] <= 0. )
	      {
		setErrmsg("ERROR: resm, powder with transverse resolution=0!") ;
		return (0) ;
	      }
	    r0 /= 2.*qsq*sqrt(rm[1][1]*rmzz) ;  
	    /* result of dky dkz integrals */
	    
	    /*
	      1/4pi K0^2 * rt2pi/sqrt(ayy) * rt2pi/sqrt(azz)
	      1/4pi K0^2 is norm for powder mosaic distribution
	      rt2pi/sqrt  are results of integrals dky dkz
	    */
	    
	    rmv = 0. ;
	  }
	else
	  {
	    den += 1. ;
	    r0 /= sqrt(1.+sv*rmzz)*sqrt(den) ;
	    rmv = .5*rmzz/(1.+sv*rmzz) ;
	    /* this puts in 1/2 in exp(-1/2 XMX) */
	  }
	
	ss = s/den ;
	if( powder )
	  {
	    for( i=0 ; i<3 ; i++ )
	      {
		rm[i][1] = 0. ; rm[1][i] = 0. ;
	      }
	    /* b[1] = 0. ; */
	  }
	/*
	  c00 -= ss*b[1]*b[1]/4. ;
	  fac = ss*b[1] ;
	*/
	for( i=0 ; i<3 ; ++i )
	  {
	    /* b[i] -= fac*rm[1][i] ; */
	    for( j=0 ; j<3 ; ++j )
	      adiag[i][j] = rm[i][j] - ss*rm[i][1]*rm[j][1] ;
	  }
	for( i=0 ; i<3 ; ++i )
	  {
	    for( j=0 ; j<3 ; ++j ) rm[i][j] = adiag[i][j] ;
	  }

	/* put in factor of two from exp(-1/2 ...) */

	/* c00 /= 2. ; */

	for( i=0 ; i<3 ; i++ )
	  {
	    /* b[i] /= 2. ; */
	    for( j=0 ; j<3 ; j++ ) {
	      rm[i][j] /= 2. ; adiag[i][j] = rm[i][j] ;
	    }
	  }

	/*
	  Now R = *r0 exp(-c00) exp( -(X-Xo)rm(X-Xo) )
	  Now worry about calculating eigenvectors and values 
	  so diagonalize rm[][]

	  With the quadratic form as Q = X rm X + b X + C
	  change variables to X' = X - Xo  to eliminate the linear terms
	  so   Xo = -1/2 rm-1 b
	  and  Q = X' rm X' + Xo rm Xo + b Xo + C
	         = X' rm X' + 1/2 b Xo + C
	  for symmetric Hermitian rm
	  D = Uadj rm U  with unitary U columns eigenvectors
	  so rm-1 = U D-1 Uadj

	  Xo  is the new origin in X-space for the ellipsoid

	  for( i=0 ; i<3 ; i++ ) indices[i] = i ;
	  ndiag = 3 ;
	  jacobi(adiag, ndiag, ndiag, indices, eigenvals, eigenvecs, &nrot )  ;

	  find any xye offset for the reduced problem and adjust constant term
	  The resulting r0 xyeoff rm are in C&N coords, 
	  and needs to be saved for
	  CS integrals
	*/


	/*
	  convert eigenvals to HW ellipsoid
	  for( i=0 ; i<neig ; i++ ) {
	  if( fabs(eigenvals[i]) <= 0. ) {
	  strcat(errmsg," !zero eigenvalue") ;
	  return (0) ;
	  }
	  eigenvals[i] = rtln2/sqrt(eigenvals[i]) ;
	  }
	*/
	/* r0 *= exp(-c00) ; */

	/************************************
	 * 3-axis Jacobian 1/(f*ki^2*kf^3*sin(2thetas))
	 * For 3-axis multiply Jacobian by f*kfsq/ki
	 * so cross-section converted to S(Q,E) in barns/meV.
	 * For 2-axis Jacobian is 1/(k^4*sin(2thetas))
	 * i.e. Jac3-axis * sqfac = Jac2-axis
	 * cross-section is delta(ki-kf) d2s/domega
	 * the delta(ki-kf) is already used in reducing the res function
	 * so cross-section is in barns
	 ************************************************/

	ktsq = (2./3.)*epeakflux/(f/2.) ; /* from max of x^4 exp(-x^2) */
	rat = kisql/ktsq ;
	fluxperk = 1.e16*rxtotalflux*rat*rat*exp(-rat)/(twopi*ki) ;
	/* Aug 2009 bug removed mref here as already applied above */
	rat = twopi*fluxperk ;
	/*
	  r0 is in (neutron/s)(sr/k) for both cases.
	  N.B. r0 is for CS->S(Q,E) in barns/sr/meV for either 2 or 3-axis
	  normally assume 2-axis CS has delta(E) for energy integral
	  the conversion to S(Q,E) only needs to take out kf/ki
	  i.e. densstates/incidentflux
	*/
	/* rat *= kf/ki ;  kf/ki handled above with Jacobian */
	r0 *= rat ;

	/* 
	   r0 has dimensions neutron/s/A-1 
	   detector efficiency set to 1.
	   mono and anal peak reflectivities set to 1.
	   rxtotalflux is given in n/A**2/s  
	   i.e. 0.04 equivalent to 4*10^14 n/cm^2/s
	   so cross section should now be in A**2 [per meV] per steradian
	   for 3-axis.
	   factor to convert barns/meV to A**2/meV
	   1 barn = 10^-24 cm**2 = 10^-24 (10^8 A)**2 = 10^-8 A**2
	   From the cross section take out a factor of 10^24
	   to convert sample volume in cc to volume in A^3
	   This makes the overall factor 10^16
	   Put this factor in fluxperk since it is the
	   same factor for *mon below.
	   With this factor in r0 the S(Q,E) which usually has
	   (Vsample/Vcell) Fsq/cell can use cc for Vsample A^3 for Vcell
	   and barn for Fsq/cell.

	   Volume/cellVolume = Ncells
	   Note Bragg scattering also has delta(k-tau)/cellVolume
	   which is dimensionless with k in A^-3 and cellV in A^3.
	   In conclusion, for 3-axis I = R0 Int( exp(-xRMx) S(Q,E) d3deltaQ dE)
	   where S(Q,E) is in barns/meV
	   C&A Wrong! For sharp dispersion surface integrated intensity
	   is NOT just r0*volume*Speak(Q,E)!
	   For Bragg peak:
	      S(Q,E)(A^2/meV) = 10^16(A^5/barn/cc) *
               (twopi^3 delta(Q3-tau)/Vcell) (Vsample/Vcell) Fsq delta(E)
	            with Vcell in A^3, Vsample in cc, Fsq in barns.

	   For Bragg rod(2D order)
	      S(Q,E)(A^2/meV) = 10^16(A^5/barn/cc) *
               (twopi^3 delta(Q2-tau)/Vcell/c*) (Vsample/Vcell) Fsq delta(E)
	            with Vcell in A^3, Vsample in cc, Fsq in barns,
		    and c* reciprocal lattice length in A^-1
		    with c* perp to 2D plane

	   For Bragg plane(1D order)
	      S(Q,E)(A^2/meV) = 10^16(A^5/barn/cc) *
               (twopi^3 delta(Q1-tau)/Vcell/(a*b*) (Vsample/Vcell) Fsq delta(E)
	            with Vcell in A^3, Vsample in cc, Fsq in barns,
		    and a* b* reciprocal lattice length in A^-1
		    with a* b* perp to 1D line, and tau along c*


	   For 1-phonon:
	      S(Q,E)(A^2/meV) = 10^16(A^5/barn/cc) *
	       (Vsample/Vcell/8pi) Fsq gamma Bose *
	       delta(E - Ephonon(q = Q - tau))

	   where Fsq = |SUMd b_d exp(iQ.d) Qunit.polariz/sqrt(M_d/Mave)|^2
	   gamma = (Mn/Mave) Dn Q^2 / |Ephonon|
	   Bose  = (1 - exp(-beta*E))^-1 signEphonon
	*/


	/********************************************************************
	 * R(Q,E) = r0*exp(-QE*RM*QE)                                      *
	 * so that
	 * I in n/sec = r0*Integral(exp() S(Q,E)in barns/meV dE d3Q)
	 * note that r0 = R0 * f * kfsq/ki  (sqfac) in n/sec/A-1
	 * crystal peak reflectivities and detector efficiency set to 1     *
	 ********************************************************************/

	/*
	  moneff = moneffperk/ki
	  so moneffperk = e-4A-1
	  Take monitor area in cm^2 = 1e16 A^2
	*/

	fhm = mos[0]*sqrt(hcgsq[0]*hcgsq[1]/
			  (hcgsq[0] + hcgsq[1] + 4.*mosq[0]))/fabs(tom) ;

	vfacM = sqrt(vcgsq[0]*vcgsq[1]/(b0msq + vcgsq[1])) ;
	mon = monArea*moneff*(twopi*twopi/rt2pi)*mref*fluxperk*ki*fhm*vfacM ;

	/*********************************************************************
	 * mono reflectivity set to 1, and monitor area for 10cm**2 so
	 * put flux back to n/cm**2/s i.e. *10^16 again
	 * note the ki factor in Imon which cancels the one in moneff
	 *********************************************************************/


	/****************************************
	 * calc volume which is                 *
	 * (pi)**(n/2)/sqrt(det(coef.matrix))
	 * i.e.Int(exp -Ax^2) = sqrt(pi/A)
	 ****************************************/

	vol = 0. ;

	if( powder ) {
	  if( rm[0][0] > 0. ) vol = sqrt(pi/rm[0][0]) ;
	} else {
	  fac = pi ;
	  if( anal ) fac *= (fac*fac*fac) ;
	  else fac *= (fac*fac) ;
	  /*
	    NB previous detrm was without mosiac corrections
	    so we need to recompute
	  */
	  detrm = 0. ;
	  for( i=0 ; i<3 ; i++ ) {
	    i1 = (i+1)%3 ;
	    i2 = (i+2)%3 ;
	    detrm += rm[0][i]*rm[1][i1]*rm[2][i2] ;
	    detrm -= rm[0][i]*rm[1][i2]*rm[2][i1] ;
	  }
	  detrm *= rmv ;
	  if( detrm > 0. ) vol = sqrt(fac/detrm) ;
	}

	if( anal ) r0 *= aref ;  /* anal is global variable */
	tot = r0 ;
	return (1) ;
}

/*
 * example monoprog and analprog
 * here using BT7 DFM
 * please check the distances !!
 */




int mxBT7monoEfocus(int N, double k)
{
  static int Nxtal = 9 ;
  static double Di = 300. ; /* source to mono dist in cm */
  static double Df = 200. ; /* mono to sample dist in cm */
  static double Dxs = 2.25   ; /* crystal spacing in cm */
  static double Dxw = 2. ;  /* xtal blade width for info */

  static double dsp = 3.35416 ; /* PG002 dsp */

  double twothetaM, sintheta, arrayAngle, rx ;

  if( N < 1 ) return Nxtal ;
  if( N > Nxtal ) N = Nxtal ;

  /* set the flags for offset xtal calc, position mode */
  offseton = 1 ;
  monomode = 1 ;

  sintheta = PI/dsp/k ;
  if( sintheta > 1. ) return 0 ;
  twothetaM = 2.*asin(sintheta) ;
  /* for energy focus */
  arrayAngle = atan(sin(twothetaM)/(Df/Di + cos(twothetaM))) ;
  rx = (0.5*(double)N - (double)Nxtal)*Dxs ;
  xposition[1] = rx*cos(arrayAngle) ;
  yposition[1] = rx*sin(arrayAngle) ;
  zposition[1] = 0. ;
  return N ;
}
int mxBT7analEfocus(int N, double k)
{
  static int Nxtal = 13 ;
  static double Di = 200. ; /* source to mono dist in cm */
  static double Df = 40. ; /* mono to sample dist in cm */
  static double Dxs = 2.25   ; /* crystal spacing in cm */
  static double Dxw = 2. ;  /* xtal blade width for info */

  static double dsp = 3.35416 ; /* PG002 dsp */

  double twothetaA, sintheta, arrayAngle, rx ;

  if( N < 1 ) return Nxtal ;
  if( N > Nxtal ) N = Nxtal ;

  /* set the flags for offset xtal calc, position mode */
  offseton = 1 ;
  analmode = 1 ;

  sintheta = PI/dsp/k ;
  if( sintheta > 1. ) return 0 ;
  twothetaA = 2.*asin(sintheta) ;
  /* for energy focus */
  arrayAngle = atan(sin(twothetaA)/(Df/Di + cos(twothetaA))) ;
  rx = (0.5*(double)N - (double)Nxtal)*Dxs ;
  xposition[1] = rx*cos(arrayAngle) ;
  yposition[1] = rx*sin(arrayAngle) ;
  zposition[1] = 0. ;
  return N ;
}


/*****************************************************
 * kif calculates ki and kf from efix, e and fixtype *
 *****************************************************/

static int kif( int ityp, double efixl, double e, double *ki, double *kf )
{
	double c=.4825924584 ;	/* 2massneutron/hbar**2 */
	double en, eni, enf  ;

	en = e ;
	if( ityp != 2 && ityp != 3 ) en = 0. ;
	if( ityp == 3 )
		{ eni = efixl + en ; enf = efixl ; }
	else
		{ eni = efixl ; enf = efixl - en ; }

	if( enf <= 0. || eni <= 0. )
	  {
	    setErrmsg("ERROR: res3 kif invalid neutron energy!") ;
	    return (0) ;
	  }
	*ki = sqrt(c*eni) ; *kf = sqrt(c*enf) ;
	return(1) ;
}



static double recvec( double hkl[], double uni[] )
{
  /* globals ast bst cst */
	int i ;

	for( i=0 ; i<3 ; ++i )
		uni[i] = hkl[0]*ast[i] + hkl[1]*bst[i] + hkl[2]*cst[i] ;
	return (unitvec(uni)) ;
}
/* reciprocal lattice functions */
static int recipCalc(double *lat, double *ang,
		     double *av, double *bv, double *cv,
		     double *as, double *bs, double *cs, double *cvol)
{
  int i ;
  double dang[3], prod[3], dval    ;
  
  for( i=0 ; i<3 ; ++i )
    {
      if( lat[i] <= 0. || ang[i] == 0. ) return(0) ;
      else    dang[i] = DEGTORAD * ang[i] ;
    }
  
  /* x-axis along latt[0], bvec in x-y plane */
  /* find components of real space vectors */
  
  av[0] = lat[0] ; av[1] = 0. ; av[2] = 0. ;
  bv[0] = lat[1]*cos(dang[2]) ;
  bv[1] = lat[1]*sin(dang[2]) ; bv[2] = 0. ;
  if( bv[1] == 0. ) return(0) ;
  cv[0] = lat[2]*cos(dang[1])  ;
  cv[1] = (lat[2]*lat[1]*cos(dang[0])-cv[0]*bv[0])/bv[1] ;
  dval = lat[2]*lat[2] - cv[0]*cv[0] - cv[1]*cv[1] ;
  if( dval <= 0. ) return(0) ;
  cv[2] = sqrt(dval) ;
  
  vecpro( bv, cv, prod ) ;
  dval = dotpro( av, prod );
  *cvol = dval ;
  
  for( i=0 ; i<3 ; ++i ) as[i] = TWOPI*prod[i]/dval ;
  vecpro( cv, av, prod ) ;
  for( i=0 ; i<3 ; ++i ) bs[i] = TWOPI*prod[i]/dval ;
  vecpro( av, bv, prod ) ;
  for( i=0 ; i<3 ; ++i ) cs[i] = TWOPI*prod[i]/dval ;
  return(1) ;
}

static int recip()
{
  int i ;
  
  if( ! recipCalc( latt, angl, avec, bvec, cvec, ast, bst, cst, &cellvol ) )
    return(0) ;

  dpstar[0] = dotpro(ast,ast) ;
  dpstar[1] = dotpro(ast,bst) ;
  dpstar[2] = dotpro(ast,cst) ;
  dpstar[3] = dotpro(bst,bst) ;
  dpstar[4] = dotpro(bst,cst) ;
  dpstar[5] = dotpro(cst,cst) ;
  
  rlat[0] = sqrt(dpstar[0]) ;
  rlat[1] = sqrt(dpstar[3]) ;
  rlat[2] = sqrt(dpstar[5]) ;
  for( i=0 ; i<3 ; ++i ) if(rlat[i] <= 0.) return(0) ;
  
  rang[0] = acos(dpstar[1]/(rlat[0]*rlat[1]))/DEGTORAD ;
  rang[1] = acos(dpstar[4]/(rlat[1]*rlat[2]))/DEGTORAD ;
  rang[2] = acos(dpstar[2]/(rlat[2]*rlat[0]))/DEGTORAD ;
  
  L[0][0] = dpstar[0] ;
  L[0][1] = dpstar[1] ;
  L[0][2] = dpstar[2] ;
  L[1][1] = dpstar[3] ;
  L[1][2] = dpstar[4] ;
  L[2][2] = dpstar[5] ;
  L[1][0] = L[0][1] ;
  L[2][0] = L[0][2] ;
  L[2][1] = L[1][2] ;
  
  return(1) ;
}
static int recipR(Conv3_resinfo *rp)
{
  /* NO globals version */
  int i ;
  double dang[3], prod[3], dval    ;
  static double twopi = 6.2831853 ;
  static double degtorad = 0.017453293 ;

  for( i=0 ; i<3 ; ++i )
    {
      if( rp->latt[i] <= 0. || rp->angl[i] == 0. ) return(0) ;
      else    dang[i] = degtorad * rp->angl[i] ;
    }

  /* x-axis along latt[0], bvec in x-y plane */
  /* find components of real space vectors */

  rp->avec[0] = rp->latt[0] ; rp->avec[1] = 0. ; rp->avec[2] = 0. ;
  rp->bvec[0] = rp->latt[1]*cos(dang[2]) ;
  rp->bvec[1] = rp->latt[1]*sin(dang[2]) ; rp->bvec[2] = 0. ;
  if( rp->bvec[1] == 0. ) return(0) ;
  rp->cvec[0] = rp->latt[2]*cos(dang[1])  ;
  rp->cvec[1] = (rp->latt[2]*rp->latt[1]*cos(dang[0]) -
		 rp->cvec[0]*rp->bvec[0])/rp->bvec[1] ;
  dval = rp->latt[2]*rp->latt[2] - 
    rp->cvec[0]*rp->cvec[0] - rp->cvec[1]*rp->cvec[1];
  if( dval <= 0. ) return(0) ;
  rp->cvec[2] = sqrt(dval) ;

  vecpro( rp->bvec, rp->cvec, prod ) ;
  dval = dotpro( rp->avec, prod );

  for( i=0 ; i<3 ; ++i ) rp->ast[i] = twopi*prod[i]/dval ;
  vecpro( rp->cvec, rp->avec, prod ) ;
  for( i=0 ; i<3 ; ++i ) rp->bst[i] = twopi*prod[i]/dval ;
  vecpro( rp->avec, rp->bvec, prod ) ;
  for( i=0 ; i<3 ; ++i ) rp->cst[i] = twopi*prod[i]/dval ;

  rp->L[0][0] = dotpro(rp->ast,rp->ast) ;
  rp->L[0][1] = dotpro(rp->ast,rp->bst) ;
  rp->L[0][2] = dotpro(rp->ast,rp->cst) ;
  rp->L[1][1] = dotpro(rp->bst,rp->bst) ;
  rp->L[1][2] = dotpro(rp->bst,rp->cst) ;
  rp->L[2][2] = dotpro(rp->cst,rp->cst) ;

  rp->L[1][0] = rp->L[0][1] ;
  rp->L[2][0] = rp->L[0][2] ;
  rp->L[2][1] = rp->L[1][2] ;

  return(1) ;
}


/* m11*x**2 + m22*y**2 + 2*m12*x*y = p**2 */

static int ellax2( double **m, double p, double hw[2], int ix, double ***pvec )
{

  /*
    Jul 98 RWE  factor or 1/2 in hw is removed since we are already
     solving the hw equation with xRx = ln2  so hw point is like sqrt(ln2/R)
     sqrt(ln2) is the factor passed in p
  */
	double diff, thet, c, s, c2, s2, sico ;
	static double PI4 = 0.78539906 ;

	if( (m[0][0]*m[1][1] - m[0][1]*m[0][1]) <= 0. ) /* degen to a line */
		{
		hw[0] = p/2/sqrt(m[0][0]*m[0][0] + m[0][1]*m[0][1]) ;
		hw[1] = p/2/sqrt(m[1][0]*m[1][0] + m[1][1]*m[1][1]) ;
	
		pvec[ix][0][0] = m[0][0]*p/2 ; pvec[ix][0][1] = m[0][1]*p/2 ;
		pvec[ix][1][0] = m[1][0]*p/2 ; pvec[ix][1][1] = m[1][1]*p/2 ;
		return (1) ;
		}

	diff = m[1][1] - m[0][0] ;
	if( diff == 0. )
		{
		if( m[0][1] == 0. )
			thet = 0. ;
		else
			thet = PI4 ;
		}
	else
		thet = atan2(2.*m[0][1],diff)/2 ;

	c = cos(thet) ; s = sin(thet) ; c2 = c*c ; s2 = s*s ; sico = 2*s*c ;

	hw[0] = p/sqrt(m[0][0]*c2 + m[1][1]*s2 - m[0][1]*sico) ;
	hw[1] = p/sqrt(m[0][0]*s2 + m[1][1]*c2 + m[0][1]*sico) ;

	pvec[ix][0][0] = hw[0]*c ; pvec[ix][0][1] = -hw[0]*s ;
	pvec[ix][1][0] = hw[1]*s ; pvec[ix][1][1] =  hw[1]*c ;

	return(1) ;
}

static int projidire ( double **r, double dir[2], double **s, double **p )
{

  /*
    s and p return the intersection and projection principal axes hw
    in terms of ehat and dir
    i.e. the ellipse for s is
    cos(thet)(s00 dir + s01 ehat) + sin(thet)(s10 dir + s11 ehat)
    
    calc intersection and projection for ellipsoid, r,
    with plane containing e-axis and
    dir in qx-qy plane
  */

  static double fact = 0.8325545 ;
  /* sqrt(ln2)/sqrt(D)  gets half width axes from D x*x = ln2 locus */
  int i, j ;
  double nx, ny, px2, py2, v1, v2, d ;

  double hw[2] ;

  static double **ri = (double**)0 ;
  static double ***ax = (double***)0 ;

  if( ri == (double**)0 )
    {
      ri = getmatrix(ri, 0, 2, 2) ;
      ax = getmatrix3(ax, 0, 0, 3, 2, 2) ;
    }
  px2 = dir[0]*dir[0] ;
  py2 = dir[1]*dir[1] ;
  ri[0][0] = px2*r[0][0] + 2.*dir[0]*dir[1]*r[0][1] + py2*r[1][1] ;
  ri[0][1] = dir[0]*r[0][2] + dir[1]*r[1][2] ;
  ri[1][0] = ri[0][1] ;
  ri[1][1] = r[2][2] ;
  if( ! ellax2( ri,fact,hw,0,ax ) ) {
    getmatrix(ri, 2, 0, 0) ;
    getmatrix3(ax, 3, 2, 0, 0, 0) ;
    ri = NULL ;
    ax = NULL ;
    return (0) ;
  }
  for( i=0 ; i<2 ; i++ )
    for( j=0 ; j<2 ; j++ ) { s[i][j] = ax[0][i][j] ; }

  /* for projection we need the unit normal to the plane */
  if( fabs(dir[0]) > fabs(dir[1]) ) {
    ny = 1./sqrt(1.+py2/px2) ; nx = -dir[1]*ny/dir[0] ;
  } else if( fabs(dir[1]) <= 0.) {
    getmatrix(ri, 2, 0, 0) ;
    getmatrix3(ax, 3, 2, 0, 0, 0) ;
    ri = NULL ;
    ax = NULL ;
    return (0) ;
  } else { nx = 1./sqrt(1.+px2/py2) ; ny = - dir[0]*nx/dir[1] ; }

  d = nx*nx*r[0][0] + 2.*nx*ny*r[0][1] + ny*ny*r[1][1] ;
  v1 = nx*dir[0]*r[0][0] + (nx*dir[1]+ny*dir[0])*r[0][1] + ny*dir[1]*r[1][1] ;
  v2 = nx*r[0][2] + ny*r[1][2] ;
  ri[0][0] -= v1*(v1/d) ;
  ri[0][1] -= v1*(v2/d) ;
  ri[1][0] = ri[0][1] ;
  ri[1][1] -= v2*(v2/d) ;
  if( ! ellax2( ri,fact,hw,0,ax ) ) {
    getmatrix(ri, 2, 0, 0) ;
    getmatrix3(ax, 3, 2, 0, 0, 0) ;
    ri = NULL ;
    ax = NULL ;
    return (0) ;
  }
  for( i=0 ; i<2 ; i++ )
    for( j=0 ; j<2 ; j++ ) { p[i][j] = ax[0][i][j] ; }
  getmatrix(ri, 2, 0, 0) ;
  getmatrix3(ax, 3, 2, 0, 0, 0) ;
  ri = NULL ;
  ax = NULL ;
  return (1) ;
}


/*
  calc projections and intersections of 3-D ellipsoid with planes and axes
  ellax2 returns principal hw lengths and hw vectors
*/
static int proji( double **r, double axpro[3], double ***s, double ***p )
{
	int j,k,l ;
	double ra,rc,D,Akj,Alj ;
	double hw[2] ;
	static double fact = 0.8325545 ; 
	/* sqrt(ln2)/sqrt(D) factor to get 1/2widthaxes from Dx*x = ln2 locus*/

	static double **ri = (double**)0 ;
	
	if( ri == (double**)0 ) ri = getmatrix(ri, 0, 2, 2) ;

	/* get hwhm from Vmatrixelem so 1/2 = exp(-zRz) */
	if( rmv > 0 ) vres = fact/sqrt(rmv) ;

	for( j=0 ; j<3 ; ++j )
	  if( r[j][j] <= 0. ) {
	    getmatrix(ri, 2, 0, 0) ;
	    ri = NULL ;
	    return (0) ;
	  }
	for( j=0 ; j<3 ; ++j )  /* j-axis and j-plane */
	  {
	    k = (j+1)%3 ;  l = (j+2)%3 ;
	    D = r[k][k]*r[l][l] - r[k][l]*r[l][k] ;
	    Akj = (r[k][l]*r[l][j] - r[k][j]*r[l][l])/D ;
	    Alj = (r[k][l]*r[k][j] - r[l][j]*r[k][k])/D ;
	    rc = r[j][j] + 2.*r[j][k]*Akj + 2.*r[j][l]*Alj
	      + r[k][k]*Akj*Akj + r[l][l]*Alj*Alj + 2.*r[k][l]*Akj*Alj ;
	    if( rc <= 0. ) 	axpro[j] = 0. ;
	    else			axpro[j] = fact/sqrt(rc) ;
	    
	    ri[0][0] = r[k][k] ; ri[0][1] = r[k][l] ;
	    ri[1][0] = r[l][k] ; ri[1][1] = r[l][l] ;
	    if( ! ellax2( ri,fact,hw,j,s ) ) {
	      getmatrix(ri, 2, 0, 0) ;
	      ri = NULL ;
	      return (0) ;
	    }
	    ra = r[j][k]/r[j][j] ;
	    ri[0][0] = r[k][k] - ra*r[k][j] ;
	    ri[0][1] = r[k][l] - ra*r[l][j] ; ri[1][0] = ri[0][1] ;
	    ri[1][1] = r[l][l] - (r[l][j]/r[j][j])*r[l][j] ;
	    if( ! ellax2( ri,fact,hw,j,p ) ) {
	      getmatrix(ri, 2, 0, 0) ;
	      ri = NULL ;
	      return (0) ;
	    }
	  }
	getmatrix(ri, 2, 0, 0) ;
	ri = NULL ;
	return (1) ;
}
		    


static void lorentz( double k, double *lorfac, double *rlorfac )
{
	double q, sn, sn2, cs, cs2 ;
	double x, x2, y, y2, suni[3], uni[3], yunil[3] ;

	q = recvec( hkle, uni ) ;
	if( step[3] != 0. || k<= 0. || q >= 2*k )
		{ *rlorfac = 0. ; *lorfac = 0. ; return ; }

	sn  = q/k/2  ;	sn2 = sn*sn ;
	cs2 = 1.-sn2 ;	cs  = sqrt(cs2) ;

	recvec( step, suni ) ;
	x = dotpro( suni, uni ) ;	x2 = x*x ;
	vecpro( zuni, uni, yunil ) ;
	y = dotpro( suni, yunil ) ;	y2 = y*y ;

	*rlorfac = tot/sqrt(rm[0][0]*x2 + 2*rm[0][1]*x*y + rm[1][1]*y2) ;
	
	if( ifix == 0 ) *rlorfac /= cs ;	/* angle Jacobian correction */

	/* isnt the powder cone correction 1/sin(twotheta) NOT 1/sin(theta)? */
	if( ifix == 1 )		/* powder with equal q-steps */
		*lorfac = 1./sn2 ;
	else if( ifix == 0 )	/* powder with equal twotheta steps */
		*lorfac = 1./cs/sn2 ;
	else
		*lorfac = 1./sqrt(sn2*x2 + 2*sn*cs*x*y + cs2*y2) ;

}

static void matmat33( M33 a, M33 b, M33 c )
{
  int i, j ;
  for( i=0 ; i<3 ; i++ ) {
    for( j=0 ; j<3 ; j++ ) {
      c[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j] + a[i][2]*b[2][j] ;
    }
  }
}
static void transpose( M33 a, M33 t )
{
  int i, j ;
  for( i=0 ; i<3 ; i++ ) {
    for( j=0 ; j<3 ; j++ ) {
      t[j][i] = a[i][j] ;
    }
  }
}


static int set_recip()
     /*
       complete rh coordinate system with z-unit up
       also calc apub = ahkl-perp-lower-omega-angle in spectrometer sense
     */
{
  int idt, i ;
  double twopi ;
  double small, asmall ;
  double MXT[3][3] ;
  double cXa, sXa, cXb, sXb, cD, sD, cD2, sD2, sDcD ;
  twopi = TWOPI ;
  
  /* set up sample coordinate system and recip space */
  idt = recip() ;
  if( !idt )
    {
      setErrmsg("ERROR: set_recip failed to set up reciprocal space!") ;
      return (0) ;
    }
  /*
    by definition buni is to lower angle side of auni
    so that going to bhkl requires higher angle in omega
    than auni
  */
  amag = recvec(ahkl,auni) ;  /* auni is in sample coordinates */
  bmag = recvec(bhkl,buni) ;  /* buni is in sample coordinates */
  if( spin[1] == -1 )          /* right scattering */
    crosspro(auni,buni,zuni) ;
  else                        /* left scattering */
    crosspro(buni,auni,zuni) ;
  unitvec(zuni) ;
  /* zuni up by definition */
  /* complete the rhc system with apun */
  crosspro(zuni,auni,apun) ;

  /* setup columns of M as hkl of unit rhc ahkl ahklperp z */
  M[0][0] = dotpro(avec, auni)/twopi ;
  M[1][0] = dotpro(bvec, auni)/twopi ;
  M[2][0] = dotpro(cvec, auni)/twopi ;
  M[0][1] = dotpro(avec, apun)/twopi ;
  M[1][1] = dotpro(bvec, apun)/twopi ;
  M[2][1] = dotpro(cvec, apun)/twopi ;
  M[0][2] = dotpro(avec, zuni)/twopi ;
  M[1][2] = dotpro(bvec, zuni)/twopi ;
  M[2][2] = dotpro(cvec, zuni)/twopi ;
  
  /* compute the global U matrix which is f transpose from my notes */
  cXa = cos(chiA) ; sXa = sin(chiA) ;
  cXb = cos(chiB) ; sXb = sin(chiB) ;
  cD  = cos(angleA) ; sD = sin(angleA) ;
  cD2 = cD*cD ; sD2 = sD*sD ; sDcD = sD*cD ;
  U[0][0] = cD2*cXa + sD2*cXb + sDcD*sXa*sXb ;
  U[0][1] = sDcD*(cXa - cXb) + sD2*sXa*sXb ;
  U[0][2] = cD*sXa - sD*cXa*sXb ;
  U[1][0] = sDcD*(cXa - cXb) - cD2*sXa*sXb ;
  U[1][1] = cD2*cXb + sD2*cXa - sDcD*sXa*sXb ;
  U[1][2] = sD*sXa + cD*cXa*sXb ;
  U[2][0] = sD*sXb - cD*sXa*cXb ;
  U[2][1] = -cD*sXb - sD*sXa*cXb ;
  U[2][2] = cXa*cXb ;
  
  matmat33( M, U, MX ) ;
  transpose( MX, MXT ) ;
  matmat33( MXT, L, MI ) ;

  /*
    columns of MX are new hklauni hklapuni
  */

  for( i=0 ; i<3 ; i++ ) {
    ahklx[i] = MX[i][0] ;
    hkla[i] = ahklx[i] ;
    phklx[i] = MX[i][1] ;
    hklp[i] = phklx[i] ;
    zhklx[i] = MX[i][2] ;
    hklz[i] = zhklx[i] ;
  }

  amag = recvec(ahklx,aunix) ;  /* auni is in sample coordinates */
  recvec(phklx,apunx) ;
  crosspro(aunix,apunx,zunix) ;

  if( spin[1] == -1 )          /* right scattering */
    crosspro(zunix,aunix,bunix) ;
  else                        /* left scattering */
    crosspro(aunix,zunix,bunix) ;
  unitvec(bunix) ;
  bhklx[0] = dotpro(avec,bunix)/twopi ;
  bhklx[1] = dotpro(bvec,bunix)/twopi ;
  bhklx[2] = dotpro(cvec,bunix)/twopi ;

  bmag = recvec(bhklx,bunix) ;  /* bunix is in sample coordinates */

  apbdot = dotpro(apunx,bunix) ;
  bmagp = fabs(bmag * apbdot) ;
  abdot  = dotpro(aunix,bunix) ;
  /* find perp vector to auni which requires
     higher omega angle to be set, i.e. apub */
  for( i=0 ; i<3 ; i++ ) apub[i] = apunx[i] ;
  if( spin[1] == 1 ) negvec(apub) ;
  /* setup hkl for apub */
  phkl[0] = dotpro(avec,apub)/twopi ;
  phkl[1] = dotpro(bvec,apub)/twopi ;
  phkl[2] = dotpro(cvec,apub)/twopi ;
  /* find the smallest nonzero index and scale to set it 1 */
  small = 1.e32 ;
  for( i=0 ; i<3 ; i++ ) {
    asmall = fabs(phkl[i]) ;
    if( asmall > 1.e-2 && asmall < small ) small = asmall ;
  }
  for( i=0 ; i<3 ; i++ ) {
    phkln[i] = phkl[i]/small ;
    if( fabs(phkln[i]) < 1.e-3 ) phkln[i] = 0. ;
  }
  
  /* setup hkl for apun and zuni the rhc
  hklp[0] = dotpro(avec,apun)/twopi ;
  hklp[1] = dotpro(bvec,apun)/twopi ;
  hklp[2] = dotpro(cvec,apun)/twopi ;
  hklz[0] = dotpro(avec,zuni)/twopi ;
  hklz[1] = dotpro(bvec,zuni)/twopi ;
  hklz[2] = dotpro(cvec,zuni)/twopi ;
  */
  
  /*
    currently xaxis is ahkl for r-scatt and apub for left
    ? should we change plotting to always use
    x-axis as ahkl and y-axis apub... NO!!
    recall apub is near bhkl but perp to ahkl
    when abflag is set for calc_resm
    xmag, ymag are used to divide plotting stuff by to get reduced units
  */
  
  if( spin[1] == -1 ) /* right scatt auni is x, apub is y */
    {
      for( i=0 ; i<3 ; i++ )
	{
	  xuni[i] = aunix[i] ;
	  yuni[i] = apub[i] ;
	}
      xmag = amag ;
      ymag = bmagp ;
    }
  else
    {
      for( i=0 ; i<3 ; i++ )
	{
	  xuni[i] = apub[i] ;
	  yuni[i] = aunix[i] ;
	}
      xmag = bmagp ;
      ymag = amag ;
    }
  

  return (1) ;
}

static double recvecR( Conv3_resinfo *rp, double hkl[], double uni[] ) ;

static int set_recipR(Conv3_resinfo *rp)
     /*
       NO globals version
       complete rh coordinate system with z-unit up
       also calc apub = ahkl-perp-lower-omega-angle in spectrometer sense
     */
{
  int i ;
  double MXT[3][3] ;
  double cXa, sXa, cXb, sXb, cD, sD, cD2, sD2, sDcD ;
  static double twopi = 6.2831853 ;
  
  /* set up sample coordinate system and recip space */
  if( ! recipR(rp) ) return (0) ;

  /*
    by definition buni is to lower angle side of auni
    so that going to bhkl requires higher angle in omega
    than auni
  */
  rp->amag = recvecR(rp,rp->ahkl,rp->auni) ;/* auni is in sample coordinates */
  rp->bmag = recvecR(rp,rp->bhkl,rp->buni) ;/* buni is in sample coordinates */
  if( rp->scattsideright )          /* right scattering */
    crosspro(rp->auni,rp->buni,rp->zuni) ;
  else                        /* left scattering */
    crosspro(rp->buni,rp->auni,rp->zuni) ;
  unitvec(rp->zuni) ;
  /* zuni up by definition */
  /* complete the rhc system with apun */
  crosspro(rp->zuni,rp->auni,rp->apun) ;

  /* setup columns of M as hkl of unit rhc ahkl ahklperp z */
  rp->M[0][0] = dotpro(rp->avec, rp->auni)/twopi ;
  rp->M[1][0] = dotpro(rp->bvec, rp->auni)/twopi ;
  rp->M[2][0] = dotpro(rp->cvec, rp->auni)/twopi ;
  rp->M[0][1] = dotpro(rp->avec, rp->apun)/twopi ;
  rp->M[1][1] = dotpro(rp->bvec, rp->apun)/twopi ;
  rp->M[2][1] = dotpro(rp->cvec, rp->apun)/twopi ;
  rp->M[0][2] = dotpro(rp->avec, rp->zuni)/twopi ;
  rp->M[1][2] = dotpro(rp->bvec, rp->zuni)/twopi ;
  rp->M[2][2] = dotpro(rp->cvec, rp->zuni)/twopi ;
  
  /* compute the global U matrix which is f transpose from my notes */
  cXa = cos(rp->chiA) ; sXa = sin(rp->chiA) ;
  cXb = cos(rp->chiB) ; sXb = sin(rp->chiB) ;
  cD  = cos(rp->angleA) ; sD = sin(rp->angleA) ;
  cD2 = cD*cD ; sD2 = sD*sD ; sDcD = sD*cD ;
  rp->U[0][0] = cD2*cXa + sD2*cXb + sDcD*sXa*sXb ;
  rp->U[0][1] = sDcD*(cXa - cXb) + sD2*sXa*sXb ;
  rp->U[0][2] = cD*sXa - sD*cXa*sXb ;
  rp->U[1][0] = sDcD*(cXa - cXb) - cD2*sXa*sXb ;
  rp->U[1][1] = cD2*cXb + sD2*cXa - sDcD*sXa*sXb ;
  rp->U[1][2] = sD*sXa + cD*cXa*sXb ;
  rp->U[2][0] = sD*sXb - cD*sXa*cXb ;
  rp->U[2][1] = -cD*sXb - sD*sXa*cXb ;
  rp->U[2][2] = cXa*cXb ;
  
  matmat33( rp->M, rp->U, rp->MX ) ;
  transpose( rp->MX, MXT ) ;
  matmat33( MXT, rp->L, rp->MI ) ;

  /*
    columns of MX are new hklauni hklapuni
  */

  for( i=0 ; i<3 ; i++ ) {
    rp->ahklx[i] = rp->MX[i][0] ;
    rp->phklx[i] = rp->MX[i][1] ;
  }

  recvecR(rp,rp->ahklx,rp->aunix) ;  /* auni is in sample coordinates */
  recvecR(rp,rp->phklx,rp->apub) ;

  /* find perp vector to auni which requires
     higher omega angle to be set, i.e. apub */
  if( ! rp->scattsideright ) negvec(rp->apub) ;

  return (1) ;
}


static char *errstr[] = {
" monochromator angle limits exceded!",
" monochromator angle limits exceded!",
" monochromator angle limits exceded!",
" sample omega  angle limits exceded!",
" sample twotheta angle limits exceded!",
" analyzer angle limits exceded!",
" analyzer angle limits exceded!",
" unable to close scattering triangle!",
" failed to set up reciprocal lattice!",
" invalid resolution ellipsoid!",
" sample magnitude of q limits exceded!",
" no step defined for scan!",
" too many iterations!" } ;




/*
  calc mono sample and anal angles in radians for given scan point 
  N.B. doesnt modify any globals
*/

static double recvecR( Conv3_resinfo *rp, double hkl[], double uni[] )
{
  int i ;
  
  for( i=0 ; i<3 ; ++i )
    uni[i] =
      hkl[0]*rp->ast[i] + hkl[1]*rp->bst[i] + hkl[2]*rp->cst[i] ;
  return (unitvec(uni)) ;
}

static double recvecST( double as[], double bs[], double cs[],
			double hkl[], double uni[] )
{
  /* globals ast bst cst */
	int i ;

	for( i=0 ; i<3 ; ++i )
		uni[i] = hkl[0]*as[i] + hkl[1]*bs[i] + hkl[2]*cs[i] ;
	return (unitvec(uni)) ;
}


static int sampangles(double scan[], double as[], double bs[], double cs[],
		      double au[], double ap[], double ei, double ef,
		      double ang[])
{
  double ki, kf, ki2, kf2 ;
  double qu[3], qv[3], qm, q2 ;
  double arg1, arg2, del1, alpha ;
  int i ;
  static double SMALL = 1.e-12 ;
  static double pi = 3.141592654 ;
  static double twopi = 6.283185308 ;
  static double Dn = 2.072141789 ;

  ki = sqrt(ei/Dn) ;
  kf = sqrt(ef/Dn) ;
  ki2 = ki*ki ;
  kf2 = kf*kf ;

  qm = recvecST(as,bs,cs,scan,qu) ;
  q2 = qm*qm ;
  arg1 = (ki2 + kf2 - q2)/(2.*ki*kf+SMALL) ;
  if( fabs(arg1) <= 1. )	ang[2] = acos(arg1) ;
  else			        ang[2] = 0. ;

  /* del1 = cos(angle between ki and Q) */
  arg1 = (ki2 + q2 - kf2)/(2.*ki*qm+SMALL) ;
  if( fabs(arg1) <= 1. )	del1 = acos(arg1) ;
  else
    { ang[1] = 0. ; return 0 ; }
  /* now find the angle between Q and ahkl */
  for( i=0 ; i<3 ; i++ ) qv[i] = qm * qu[i] ;
  arg1 = dotpro(qv,au) ;
  /* Jul 1998 removed spin[1] factor from arg2
     since now apub is always to lower omega angle */
  arg2 = dotpro(qv,ap) ;
  if( fabs(arg1) < 1.e-6 && fabs(arg2) < 1.e-6 )
    { alpha = 0. ; }
  else
    { alpha = atan2(arg1,arg2) ; }
  /*
    alpha is angle of Q [-pi,pi] in x=apub y=auni coordinate sys
    so subtract PIOVER2 to get angle Q makes with ahkl 
    with sign sense of omega
    alpha2 = alpha - PIOVER2 ;
    Now the angle ahkl makes wrt KI is del1 + alpha2 
    with positive towards the scattering side.
    the zero position for ahkl is at +PIOVER2 from KI, 
    again with pos toward the scat side.
    so the angle ahkl makes wrt its zero position in the same sense 
    is del1 + alpha2 - PIOVER2.
    But this is the opposite sense of omega, so take the negative to give,
    PIOVER2 - del1 - alpha2 = PIOVER2 - del1 - alpha + PIOVER2 
    = PI - del1 - alpha
  */
  ang[1] = fmod(pi-del1-alpha,twopi) ;
  while( ang[1] >= twopi ) ang[1] -= twopi ;
  while( ang[1] < 0. )     ang[1] += twopi ;
  return 1 ;
}
static int sampanglesK(double scan[], double as[], double bs[], double cs[],
		      double au[], double ap[], double ki, double kf,
		      double ang[])
{
  double ki2, kf2 ;
  double qu[3], qv[3], qm, q2 ;
  double arg1, arg2, del1, alpha ;
  int i ;
  static double SMALL = 1.e-12 ;
  static double pi = 3.141592654 ;
  static double twopi = 6.283185308 ;

  ki2 = ki*ki ;
  kf2 = kf*kf ;

  qm = recvecST(as,bs,cs,scan,qu) ;
  q2 = qm*qm ;
  arg1 = (ki2 + kf2 - q2)/(2.*ki*kf+SMALL) ;
  if( fabs(arg1) <= 1. )	ang[2] = acos(arg1) ;
  else			        ang[2] = 0. ;

  /* del1 = cos(angle between ki and Q) */
  arg1 = (ki2 + q2 - kf2)/(2.*ki*qm+SMALL) ;
  if( fabs(arg1) <= 1. )	del1 = acos(arg1) ;
  else
    { ang[1] = 0. ; return 0 ; }
  /* now find the angle between Q and ahkl */
  for( i=0 ; i<3 ; i++ ) qv[i] = qm * qu[i] ;
  arg1 = dotpro(qv,au) ;
  /* Jul 1998 removed spin[1] factor from arg2
     since now apub is always to lower omega angle */
  arg2 = dotpro(qv,ap) ;
  if( fabs(arg1) < 1.e-6 && fabs(arg2) < 1.e-6 )
    { alpha = 0. ; }
  else
    { alpha = atan2(arg1,arg2) ; }
  /*
    alpha is angle of Q [-pi,pi] in x=apub y=auni coordinate sys
    so subtract PIOVER2 to get angle Q makes with ahkl 
    with sign sense of omega
    alpha2 = alpha - PIOVER2 ;
    Now the angle ahkl makes wrt KI is del1 + alpha2 
    with positive towards the scattering side.
    the zero position for ahkl is at +PIOVER2 from KI, 
    again with pos toward the scat side.
    so the angle ahkl makes wrt its zero position in the same sense 
    is del1 + alpha2 - PIOVER2.
    But this is the opposite sense of omega, so take the negative to give,
    PIOVER2 - del1 - alpha2 = PIOVER2 - del1 - alpha + PIOVER2 
    = PI - del1 - alpha
  */
  ang[1] = fmod(pi-del1-alpha,twopi) ;
  while( ang[1] >= twopi ) ang[1] -= twopi ;
  while( ang[1] < 0. )     ang[1] += twopi ;
  return 1 ;
}
static void specangles( double scan[], double ang[] )
{
  /* uses global ast, bst, cst, auni, apub, dsps[], ifix, efixd */
  double ki, kf ;
  double arg1, arg2 ;

  if( ! kif(ifix,efixd,scan[3],&ki,&kf) ) { scan[3] = 0. ; }

  arg1 = PI/(ki*dsps[0]) ;
  arg2 = PI/(kf*dsps[1]) ;
  if( arg1 > 1. )
    { ang[0] = 0. ; }
  else
    { ang[0] = 2.*asin(arg1) ; } /* in radians */
  if( arg2 > 1. )
    { ang[3] = 0. ; }
  else
    { ang[3] = 2.*asin(arg2) ; } /* in radians */

  sampanglesK(scan, ast, bst, cst, aunix, apub, ki, kf, ang) ;
}

/* calc spec angles in degrees */
static void specanglesd( double scan[], double ang[] )
{
  int i ;
  specangles(scan, ang) ;
  for( i=0 ; i<4 ; i++ ) ang[i] *= RADTODEG ;
}

static double etok( double e )
{
	return ( sqrt((fabs(e)+1.e-16)/NEUTSTIFF) ) ;
}



double geocsmin, root[2], xmn ;
int winit ;


static double countrates( double scan[] )
{
  double rate, countratelimit, sinthet, sintwothet, gmin ;
  double fac, rab, rex ;
  double ang[4] ;

  relint = 0. ;
  if( ! geocshkle( scan, &gmin ) ) return 0. ;
  /* geocshkle calcs globals: countrate and Fsqeff */
  specangles(scan, ang) ;
  sinthet = sin(ang[2]/2.) ;
  if( sinthet <= 0. ) sinthet = 1. ;
  rate = sampleArea*sampleThick*mref*aref ;


  /*
    now calc the absorption and extinction corrections
    note extinctLength is calculated in cm
  */

  /* relint calc is in geocshkle */

  rate *= Fsq ;
  if( Fsqeff <= 0. ) Fsqeff = 1.e-4 ;

  if( mosa[1] <= 0. )  /* primary extinction */
    {
      rate *= exp(-absorbCoef*sampleThick) ;
      extinctLength = cellvol*sinthet/(lamb*sqrt(Fsqeff)*10000.) ;
      fac = sampleThick/extinctLength ;
      rate *= tanh(fac)/fac ;
    }
  else /* else secondary extinction */
    {
      fac = (180.*60.)/(sqrt(TWOPI)*PI*mosa[1]) ;
      /* fac = W(0) */
      sintwothet = sin(ang[2]) ;
      extinctLength = cellvol*cellvol*sintwothet/(lamb*lamb*lamb*Fsqeff);
      /* sets extinctLength to 1/Q */
      extinctLength /= fac ;
      if(absorbCoef <= 0.)
	{
	  extinctLength *= sinthet ;
	  rate /= 1. + relint*sampleThick/extinctLength ;
	}
      else
	{
	  rex = relint/(absorbCoef*extinctLength) ;
	  rab = absorbCoef*sampleThick/sinthet ;
	  rate /= rab ;
	  fac = sqrt(1. + 2.*rex) ;
	  rate /= 1. + rex + fac/tanh(rab*fac) ;
	}
    }


  rate *= countrate ;
  if( monArea > 0. && moneff > 0. )
    {
      countratelimit = mon*(sampleArea/monArea)/moneff ;
      if( rate > countratelimit ) rate = countratelimit ;
    }
  return (rate) ;
}


/*
  calc width for simple geometric cross-sections
*/  
  
static int width_calc()
{
  double xwid, hhhhmin ;
  double xl1, xl2, x1, x2, xl, xr, xmin, xphh, xmhh, xint ;
  double vect[3], vectuni[3] ;
  double scanm[4], scanp[4], anglm[4], anglp[4], cotm, cotp ;
  
  
  double stepmax, stepmax2 ;
  int iwcsave, i, foundpeak, indmax ;
  
  /* first find max intensity for the given scan */
  
  for( i=0 ; i<4 ; i++ )
    {
      peak[i] = 0. ;
      fwhm[i] = 0. ;
      angls[i] = 0. ;
      fwhmang[i] = 0. ;
    }
  peakmag = 0. ;
  fwhmmag = 0. ;
  cota = 0. ;
  fwhmcot = 0. ;
  
  xl1 = (1. - (double)npts)/2. ;
  xl2 = -xl1 ;
  
  iwcsave = iwctyp ;	iwctyp = 1 ;	/* force dynamic for first min */
  /* improved search RWE oct 99 */
  /* estimate the min from cspoint */
  /* RWE May 2007 except for linear_dispersion use scan center as min guess */
  indmax = 0 ;
  for( i=1 ; i<5 ; i++ )
    { if( (stepmax = fabs(step[i-1])) > 0. ) { indmax = i ; break ; } }
  if( !indmax ) return (0) ;

  indmax-- ;
  for( i=indmax+1 ; i<4 ; i++ )
    { 
      if( (stepmax2 = fabs(step[i])) > stepmax )
	{ indmax = i ; stepmax = stepmax2 ; }
    }

  xmin = 0. ;
  if( icstyp < 4 && icstyp >= 0 ) {
    xmin = (point[indmax] - hkle[indmax])/step[indmax] ;
    if( xmin + 5. > xl2 ) xl2 = xmin + 6. ;
    if( xmin - 5. < xl1 ) xl1 = xmin - 6. ;
  }

  foundpeak = 0 ;
  if( icstyp < 5 && icstyp >= 0 )
    foundpeak = findmin(xl1, xl2, 1., &xmin, &x1, &x2, &geocsmin, -1., geocs);
  //printf("width_calc: foundpeak = %d xmin = %g\n", foundpeak, xmin) ;
  if( !foundpeak )
    {
      calc_resm(dopt, 0) ;
      specanglesd(dopt, angls) ;
      return 0 ;
    }
  for( i=0 ; i<4 ; i++ ) {
    peak[i]  = hkle[i] + xmin*step[i] ;
    efixd = efix + xmin*defix ;
  }
  specanglesd(peak, angls) ;
  
  xmn = 0. ;
  xwid = 0. ;
  iwctyp = 0 ; /* turn off dynamic */
  
  /* try to refine the peak position */
  if( icstyp < 5 )
    {
      /* removed mnbrak as this is replaced by bisectsearch above */
      /* if( mnbrak( &x1, &x2, &x3, &f1, &f2, &f3, geocs ) ) */
      /* countrate will be peakcountrate after this call to brent */

      /*
	Sep 2002 try replace brent with findmin with precision
      */
      /* if( ! brent( x1, xmin, x2, geocs, tol, maxit, &xmn, &geocsmin ) ) */
      xint = fabs(x1 - x2)/10. ;
      xl = x1 ;
      xr = x2 ;
      if(! findmin(xl, xr, xint, &xmin, &x1, &x2, &geocsmin, fmprec, geocs))
	{
	  setErrmsg("ERROR: failed res width solution!") ;
	  geocshkle( hkle, &geocsmin ) ;
	  return 0 ;
	}
    }
  xmn = xmin ;  
  for( i=0 ; i<4 ; i++ ) {
    peak[i]  = hkle[i] + xmn*step[i] ;
    efixd = efix + xmn*defix ;
  }
  specanglesd(peak, angls) ;
  
  for( i=0 ; i<3 ; i++ ) { vect[i] = peak[i] ; }
  peakmag = recvec(vect,vectuni) ;
  
  /* find bracket for halfheight zeros */

  xphh = xmn + 1. ;
  if( xphh > xl2 ) xphh = (xmn + xl2)/2. ;
  
  /*
    changed to bisectsearch RWE Oct 99
    if( mnbrak( &xstar, &xphh, &x3, &f1, &f2, &f3, hhhh ) )
    note stepsize decr to 0.1 for hh
  */
  if( findmin(xmn, xl2, 0.1, &xphh, &x1, &x2, &hhhhmin, -1., hhhh) )
    {
      /*
	change to brent on hhhh 
	instead of zbrent root finder RWE Oct 99
      */
      
      /* if( zbrent( hhgeocs, xstar, x3, tol, maxit, &root[1] ) ) */
      
      /* if( brent( x1, xphh, x2, hhhh, tol, maxit, &xhhp, &hhhhmin ) ) */
      xl = x1 ;
      xr = x2 ;
      xint = fabs(xr - xl)/10. ;
      if( findmin(xl, xr, xint, &xphh, &x1, &x2, &hhhhmin, fmprec, hhhh) )
	{
	  xmhh = xmn - (xphh - xmn) ;
	  if( xmhh < xl1 ) xmhh = (xmn + xl1)/2. ;
	  /* if( mnbrak( &xstar, &xmhh, &x3, &f1, &f2, &f3, hhhh ) ) */
	  if( findmin(xl1, xmn, 0.1, &xmhh, &x1, &x2, &hhhhmin, -1., hhhh) )
	    {
	      /*
		if( zbrent( hhgeocs, xstar, x3, tol, maxit, &root[0] ) )
	      */
	      /* if(brent( x1, xmhh, x2, hhhh, tol, maxit, &xhhm, &hhhhmin))*/
	      xl = x1 ;
	      xr = x2 ;
	      xint = fabs(xr - xl)/10. ;
	      if(findmin(xl, xr, xint, &xmhh, &x1, &x2, &hhhhmin,fmprec, hhhh))
		{
		  xwid = fabs(xphh - xmhh) ;
		  for( i=0 ; i<4 ; i++ )
		    {
		      fwhm[i] = xwid*step[i] ;
		      scanm[i] = peak[i] - fwhm[i]/2. ;
		      scanp[i] = peak[i] + fwhm[i]/2. ;
		    }
		  for( i=0 ; i<3 ; i++ ) { vect[i] = fwhm[i] ; }
		  fwhmmag = recvec(vect,vectuni) ;
		  efixd = efix + xmhh*defix ;
		  specanglesd(scanm, anglm) ;
		  efixd = efix + xphh*defix ;
		  specanglesd(scanp, anglp) ;
		  for( i=0 ; i<4 ; i++ ) {
		    fwhmang[i] = anglp[i] - anglm[i] ;
		  }
		  if ( anglp[3] > 0. && anglm[3] > 0. ) {
		    cotp = 1./tan( DEGTORAD * anglp[3]/2. ) ;
		    cotm = 1./tan( DEGTORAD * anglm[3]/2. ) ;
		    fwhmcot = cotp - cotm ;
		  }
		  iwctyp = iwcsave ;
		  return (1) ;
		}
	    }
	}
    }
  iwctyp = iwcsave ;
  return (0) ;
}


static int hhhh( double x, double *v )
{
	static double hh ;
	if( ! hhgeocs(x, &hh) ) return (0) ;
	*v = hh*hh ;
	return (1) ;
}

static int hhgeocs( double x, double *v )
{
        static double hh ;
	if( ! geocs(x, &hh) ) return (0) ;
        *v = 0.5*geocsmin - hh ;
	return (1) ;
}


/* calc simple geometric convolutions used by widths */
/* note that negative of intensity is returned in *v */

static int geocs( double x, double *v )
{
  int i ;
  double scan[4] ;
  /* calculate the scan vector  vec = vecs0 + x*vecsd */
  for( i=0 ; i<4 ; i++ )
    {
      scan[i]  = hkle[i]  + x*step[i]  ;   /* hkl coords   */
      efixd = efix + x*defix ;
    }
  return (geocshkle(scan, v)) ;
}

/*
  geocshkle  does convolutions of resol func 
  with simple geometric cross sections.
  x in geocs is parameterized point in scan, i.e.
  e = e0 + x*dele  etc.
  returns convolution of resolution with a simple
  geometrical cross section which can be a point, line or plane.
  There is a global flag, icstyp,
  which determines which cross section type.
  
  convert scan and cross section vectors to ahklunit and perp rhc
  Then exponent in resolution function is quadratic in the scan
  parameter x.
  xscanvec = xscan0vec + x*xscanstepvec
  cstype 0 = Bragg point   xvec = xcsvec0
  cstype 1 = Bragg line	   xvec = xcsvec0 + l*xcsdirvec1
  cstype 2 = Bragg plane   xvec = xcsvec0 + l*xcsdirvec1 + m*xcsdirvec2
  At a fixed scan point convoluting the res function with these
  simple cross sections is straight forward.

  Also calcs countrate

	   For Bragg peak:
	      S(Q,E)(A^2/meV) = 10^16(A^5/barn/cc) *
               (twopi^3 delta(Q3-tau)/Vcell) (Vsample/Vcell) Fsq delta(E)
	            with Vcell in A^3, Vsample in cc, Fsq in barns.

	   For Bragg rod(2D order)
	      S(Q,E)(A^2/meV) = 10^16(A^5/barn/cc) *
               (twopi^3 delta(Q2-tau)/Vcell/c*) (Vsample/Vcell) Fsq delta(E)
	            with Vcell in A^3, Vsample in cc, Fsq in barns,
		    and c* reciprocal lattice length in A^-1
		    with c* perp to 2D plane

	   For Bragg plane(1D order)
	      S(Q,E)(A^2/meV) = 10^16(A^5/barn/cc) *
               (twopi^3 delta(Q1-tau)/Vcell/(a*b*) (Vsample/Vcell) Fsq delta(E)
	            with Vcell in A^3, Vsample in cc, Fsq in barns,
		    and a* b* reciprocal lattice length in A^-1
		    with a* b* perp to 1D line, and tau along c*


	   For 1-phonon:
	      S(Q,E)(A^2/meV) = 10^16(A^5/barn/cc) *
	       (Vsample/Vcell/8pi) Fsq gamma Bose *
	       delta(E - Ephonon(q = Q - tau))

	   where Fsq = |SUMd b_d exp(iQ.d) Qunit.polariz/sqrt(M_d/Mave)|^2
	   gamma = (Mn/Mave) Dn Q^2 / |Ephonon|
	   Bose  = (1 - exp(-beta*E))^-1 signEphonon

*/


static int geocshkle( double scan[], double *v )
{
  int i, j, ip, jp, retval ;
  double a0, dmag, cs1mag, cs2mag, cs3mag ;
  double xdelt[4], csd1[4], csd2[4], csd3[4], delt[4] ;
  double cs3u[3], cs2u[3], cs1u[3], duni[3] ;
  double dvec[3] ;
  
  double d11, d12, d22, d1d, d2d ;
  double cphonon[3], c0phonon, Ephonon ;
  double xvec[3], det, sin12 ;
  double mvec[3], temp ;
  /*
    static double f  = 4.144283578
    hbar**2/massneutron (mev A**2)
  */
  static double Dn = 2.072141789;
  static double pi = PI ;
  
  static double **mmat = (double**)0 ;
  static double **minv = (double**)0 ;
  static double **dmat = (double**)0 ;
  static double **dinv = (double**)0 ;
  static double **rmwl = (double**)0 ;  /* local copy of rmw */
  static double **rmph = (double**)0 ;  /* for phonon calc */
  static double **Sm = (double**)0 ;
  
  static int elimpref[4] = {2, 0, 1, 3} ;
  int vsolved[4] ;
  static int vremain[4] = {1, 1, 1, 1} ;
  /* init flag all variables active */
  int nequdone, llinear ;
  
  static double rmb[4] = {0., 0., 0., 0.} ;
  double rmc ;
  
  
  static EQUATION *gequp = (EQUATION*)0 ;
  
  
  if( gequp == (EQUATION*)0 )
    {
      gequp = (EQUATION*)malloc(3*sizeof(EQUATION)) ;
      gequp[0].coef = (double*)malloc(5*sizeof(double)) ;
      gequp[1].coef = (double*)malloc(5*sizeof(double)) ;
      gequp[2].coef = (double*)malloc(5*sizeof(double)) ;
      rmwl = getmatrix(rmwl, 0, 3, 3) ;
      rmph = getmatrix(rmph, 0, 4, 4) ;
      dmat = getmatrix(dmat, 0, 3, 3) ;
      dinv = getmatrix(dinv, 0, 3, 3) ;
      mmat = getmatrix(mmat, 0, 3, 3) ;
      minv = getmatrix(minv, 0, 3, 3) ;
      Sm = getmatrix(Sm, 0, 2, 2) ;
    }
  
  
  
  /* calculate the scan vector  vec = vecs0 + x*vecsd */
  for( i=0 ; i<4 ; i++ )
    {
      delt[i]  = point[i] - scan[i]    ;   /* diff w point */
    }
  /* calc resolution matrix, in ahkl coordinates used by cs-calc */
  if( iwctyp || winit )	/* dynamic calc get new res matrix */
    {
      if( ! calc_resm( scan, 1 ) ) 	{
	*v = 0. ;
	retval = 0 ;
	goto ctrRet ;
      } else				{
	winit = 0 ; }
    }
  for( i=0 ; i<3 ; i++ )
    {
      for( j=0 ; j<3 ; j++ ) {
	rmq[i][j] = rmw[i][j] ; rmwl[i][j] = rmw[i][j] ;
	rmph[i][j] = rmw[i][j] ;
      }
      rmq[2][i] = 0. ; rmq[i][2] = 0. ;
      rmph[3][i] = 0. ; rmph[i][3] = 0. ;
      rmb[i] = 0. ;
    }
  rmq[2][2] = rmv ;
  rmph[3][3] = rmv ;
  rmb[3] = 0. ;
  
  /*
    if not dynamic rm previously calc with ahkl coords
    convert all points and directions to auni and apun coords
    NB any offset of the resolution function is returned from
    calc_resm in ahkl ahklperp-rhc
    xdelt is to be distance from resolution ellipsoid center to CSpoint
    but ellipsoid center is at xyeoffw added to scan in auni, apun coords
  */
  
  qmag = recvec(scan,quni) ;
  dmag = recvec(delt,duni) ;
  xdelt[0] = dotpro(auni,duni)*dmag - xyeoffw[0] ;
  xdelt[1] = dotpro(apun,duni)*dmag - xyeoffw[1] ;
  xdelt[2] = dotpro(zuni,duni)*dmag - qzoff ;
  xdelt[3] = delt[3] - xyeoffw[2] ;
  
  /*
    for crystal of 1cm**3 volume and structure factor of 1e-24cm
    i.e. in barns
    for point CS assume Bragg peak which is
    F^2 Vsample (2pi)^3 / cellvol^2 delta(Q3-tau3)  OR
    (|Fcell/ncell|^2) Vsample(cm^3) (2pi)^3/(cellvol/ncell)^2
    
    for Bragg rod assume independent 2D scattering planes
    so modify Bragg peak by  delta(Q2-tau2)/(twopi/zcell)
    where twopi/zcell is the recip cell dimension
    along the rod direction.
    This 1/cs1mag will cancel cs1mag from Jacobian of integral
    
    for Bragg plane assume independent 1D scattering lines
    so modify Bragg peak by  delta(Q1-tau1)/(twopi/acell * twopi/bcell)
    i.e. divide by cell area perp to 1D order.
    This 1/(cs1magxcs2mag) cancels the Jacobian except for angle factor
    
    For phonon 1cc sample      S(Q,E) = (1/vcell) Fsq gamma Bose * 1A^3
    * delta(E - Ephonon(Q))
    where Fsq = |SUMd b_d exp(iQ.d) Qunit.polariz/sqrt(M_d/Mave)|^2
    gamma = (Mn/Mave) Dn Q^2 / |Ephonon|
    Bose  = (1 - exp(-beta*E))^-1 signEphonon
    We take Bose = 1 and Fsq = 1 and Mn/Mave = 1
    typically Fsq is larger and Mn/Mave is smaller (Mn = neutron mass)
    
    General form is r0 * Ccs * pi^(n/2) exp(-Vo M Vo) |S|^-1/2 exp(X Q X)
    For Bragg pt,line,plane
    Ccs = F^2 Vsample (2pi)^3 / cellvol^2 / Qcell0,length,area
    Vo = Pcs - QEo
    S  = D-adj M D
    Q  = S^-1 adj
    X-adj = Vo M D = (Voq Mqq + Voe Me-adj) D 
    M = resol matrix
    Dij = Dcsq ithcomp, jthDcs  3xn matrix
    
    For phonon n=3
    Ccs = (Vsample/cellvol) Fsq gamma Bose
    Fsq = |Sum d  b_d exp(i Q.d) Qhat.polariz/sqrt(M_d/Mave)|^2
    gamma = (Mneut/Mave) Dn Qphonon^2 / |Ephonon|
    Bose = (1 - exp(-beta Ephonon))^-1 signEphonon
    Qphonon = Q0 = qmag
    
    Vo = (0,0,co - Eo,0)
    S = Mqq + Mee c c-adj + c Me-adj + Me c-adj
    Me_i = Mei  3-vector
    Q = S^-1 adj
    X-adj = (co - Eo) (Mee c-adj + Me-adj)
    c-adj = De-adj D^-1
    co = Pe - c-adj Pq
  */
  
  
  /* for Bragg scattering */
  
  if( icstyp < 3 )
    {
      countrate = TWOPI*TWOPI*TWOPI/(cellvol*cellvol) ;
      a0 = -vrv(xdelt,xdelt) ; /* all Bragg CS have this factor */
      
      /*
	if there is a linear constraint on Resolution
	convert it to auni apun zuni and then xdelt must satisfy
	the linear constraint 
	for a Bragg peak to lie on the 2D resol plane
	i.e. xdelt on constraint plane
      */
      
      if( icstyp == 0 ) /* Bragg point */
	{
	  if( a0 > 26. ) a0 = 26. ;
	  if( a0 < -26. ) a0 = -26. ;
	  relint = exp(a0) ;
	  countrate *= tot*relint ;
	  Fsqeff = Fsq ;
	  *v = -countrate ;
	  retval = 1 ;
	  goto ctrRet ;
	}
      
      countrate *= pow(pi,(double)icstyp/2.) ;  
      /* from integrals below */
      
      cs1mag = recvec(dire1,cs1u) ;
      csd1[0] = dotpro(auni,cs1u) ;
      csd1[1] = dotpro(apun,cs1u) ;
      csd1[2] = dotpro(zuni,cs1u) ;
      csd1[3] = 0. ;
      d11 = vrvq(csd1,csd1) ;
      d1d = vrv(csd1,xdelt) ;
      
      if( icstyp == 1 ) /* Bragg line with recip cell dimension cs1mag */
	{
	  a0 += d1d*d1d/d11 ;
	  temp = sqrt(d11) ;
	  if( a0 > 26. ) a0 = 26. ;
	  if( a0 < -26. ) a0 = -26. ;
	  relint = exp(a0) ;
	  countrate *= tot*relint/temp ;
	  Fsqeff = Fsq/temp ;
	  *v = -countrate ;
	  retval = 1 ;
	  goto ctrRet ;
	}	    
      
      /* Bragg (elastic) plane with cell area cs1mag*cs2mag*csd1xcsd2 */
      cs2mag = recvec(dire2,cs2u) ;
      csd2[0] = dotpro(auni,cs2u) ;
      csd2[1] = dotpro(apun,cs2u) ;
      csd2[2] = dotpro(zuni,cs2u) ;
      csd2[3] = 0. ;
      /* make sure the second direction is perpendicular to first */
      crosspro(csd1,csd2,csd3) ;
      Fsqeff = 0. ;
      if( (sin12 = unitvec(csd3)) <= 0. ) {
	retval = 0 ;
	goto ctrRet ;
      }
      crosspro(csd3,csd1,csd2) ;
      
      d12 = vrvq(csd1,csd2) ;
      d22 = vrvq(csd2,csd2) ;
      d2d = vrv(csd2,xdelt) ;
      dmag = d11*d22 - d12*d12 ;
      
      a0 += (d11*d2d*d2d + d22*d1d*d1d - 2.*d12*d1d*d2d)/dmag ;
      temp = sin12*sqrt(dmag) ;
      
      if( a0 > 26. ) a0 = 26. ;
      if( a0 < -26. ) a0 = -26. ;
      relint = exp(a0) ;
      countrate *= tot*relint/temp ;
      Fsqeff = Fsq/temp ;
      *v = -countrate ;
      retval = 1 ;
      goto ctrRet ;
    }
  
  if( icstyp == 3 ) /* incoherent in Q, but delta(E-Ep) */
    {
      dmag = rm[0][0]*rm[1][1] - rm[0][1]*rm[0][1] ;
      temp = rm[0][0]*rm[1][2]*rm[1][2]
	+ rm[1][1]*rm[0][2]*rm[0][2] ;
      temp -= 2.*rm[0][1]*rm[0][2]*rm[1][2] ;
      temp = rm[2][2] - temp/dmag ;
      a0 = -temp*xdelt[3]*xdelt[3] ;
      countrate = pow(TWOPI,1.5) ;
      
      if( a0 > 26. ) a0 = 26. ;
      if( a0 < -26. ) a0 = -26. ;
      relint = exp(a0) ;
      /* with rtwpi/sqrt(dmag) this is integral over X and Y */
      countrate *= tot*relint/sqrt(rmv*dmag) ;
      Fsqeff = Fsq ;
      *v = -countrate ;
      retval = 1 ;
      goto ctrRet ;
    }
  
  llinear = 0 ;
  if( icstyp == 4 || icstyp == 5 ) /* linear dispersion OR uniform */
    {
      if( icstyp == 4 )
	{
	  
	  cs1mag = recvec(dire1,cs1u) ;
	  dmat[0][0] = dotpro(auni,cs1u)*cs1mag ;
	  dmat[1][0] = dotpro(apun,cs1u)*cs1mag ;
	  dmat[2][0] = dotpro(zuni,cs1u)*cs1mag ;
	  dvec[0] = dire1[3] ;
	  
	  cs2mag = recvec(dire2,cs2u) ;
	  dmat[0][1] = dotpro(auni,cs2u)*cs2mag ;
	  dmat[1][1] = dotpro(apun,cs2u)*cs2mag ;
	  dmat[2][1] = dotpro(zuni,cs2u)*cs2mag ;
	  dvec[1] = dire2[3] ;
	  
	  cs3mag = recvec(dire3,cs3u) ;
	  dmat[0][2] = dotpro(auni,cs3u)*cs3mag ;
	  dmat[1][2] = dotpro(apun,cs3u)*cs3mag ;
	  dmat[2][2] = dotpro(zuni,cs3u)*cs3mag ;
	  dvec[2] = dire3[3] ;
	  
	  /*
	    Ecs = c0 + cvec*q(a,ap,z)
	    so dmat must be inverted to find c0 and cvec
	  */
	  
	  if( ! inv33( dmat, dinv ) ) {
	    retval = 0 ;
	    goto ctrRet ;
	  }
	  vecmat3( dvec, dinv, cphonon ) ;
	  c0phonon = point[3] - dotpro( cphonon, xdelt ) ;
	  
	  /* OK since xdelt was corrected for offset */
	  /* and recall delt = point - scan */
	  
	  Ephonon = fabs(c0phonon) ;
	  if( Ephonon < 1.e-4 ) Ephonon = 1.e-4 ;
	  /* pi^3/2 = 5.568   factors of pi done below */
	  /* this is energy of phonon at dQ = 0 */
	  
	  countrate = tot*Dn*qmag*qmag/fabs(Ephonon)/cellvol ;
	  /*
	    this is part of countrate specific to phonon
	    i.e. Dn Q^2 / Ephonon / cellvol
	    where we have taken the Bose factor * FsqPhonon/4pi = 1
	    and again for 1cc sample as above
	  */
	  
	  /*  using manual reduction from delta(Ephonon - ...)
	      edelt = c0phonon - scan[3] ;
	      a0 = -rmw[2][2] ;
	      for( i=0 ; i<2 ; i++ ) mvec[i] = rmw[2][i] ;
	      mvec[2] = 0. ;
	      for( i=0 ; i<3 ; i++ )
	      for( j=0 ; j<3 ; j++ )
	      mmat[i][j] = rmw[2][2]*cphonon[i]*cphonon[j]
	      + mvec[i]*cphonon[j] + mvec[j]*cphonon[i] ;
	      for( i=0 ; i<2 ; i++ )
	      for( j=0 ; j<2 ; j++ )
	      mmat[i][j] += rmw[i][j] ;
	      mmat[2][2] += rmv ;
	      det = det33(mmat) ;
	      if( det <= 0. ) return (0) ;
	      inv33( mmat, minv ) ;
	      for( i=0 ; i<3 ; i++ )
	      {
	      mvec[i] += rmw[2][2]*cphonon[i] ;
	      for( j=0 ; j<3 ; j++ ) mmat[i][j] = minv[j][i] ;
	      }
	      matvec3(mmat,mvec,xvec) ;
	      a0 += dotpro(mvec,xvec) ;
	      a0 *= edelt*edelt ;
	  */
	  /*
	    another technique is to use lreduce 
	    to apply the phonon delta function
	    This has the advantage that any additional 
	    linear constraint is handled in
	    the same way
	    So generate the phonon equation in gequp[linear]
	  */
	  llinear = 0 ;
	  gequp[llinear].coef[0] = scan[3] + xyeoff[2] - c0phonon ;
	  gequp[llinear].coef[1] = -cphonon[0] ;
	  gequp[llinear].coef[2] = -cphonon[1] ;
	  gequp[llinear].coef[3] = 1. ;
	  /* dE coef for delta(dE ... ) */
	  gequp[llinear].coef[4] = -cphonon[2] ;  /* z coef */
	  ++llinear ;
	  
	  gequp[0].nvar = 4 ;
	  
	  /* will eliminate elimpref[0] first etc. */
	  /* corresponding coef index = *elimpref + 1 */
	  elimpref[0] = 2 ; elimpref[1] = 0 ;
	  elimpref[2] = 1 ; elimpref[3] = 3 ;
	  vremain[0] = 1 ;  vremain[1] = 1 ;
	  vremain[2] = 1 ;  vremain[3] = 1 ; /* init all vars activ */
	}
      else  /* icstype == 5 uniform CS  with no linear equs */
	{
	  /* part of CS specific to uniform CS */
	  countrate = tot ;
	}
      if( llinear )
	{
	  /* phonon will have linear >= 1 */
	  nequdone = lreduce( llinear, gequp, elimpref,
			      vsolved, vremain, eps ) ;
	  if( nequdone < 0 )
	    {
	      *v = 0. ;
	      relint = 0. ;
	      countrate = 0. ;
	      Fsqeff = 0. ;
	      retval = 0 ;
	      goto ctrRet ;
	    }
	  
	  /* now reduce the XYEZ integrand quadratic form in rmph */
	  
	  rmc = 0. ;
	  qreducem( rmph, rmb, &rmc, 4, llinear,
		    vsolved, vremain, gequp ) ;
	  a0 = -rmc ;
	  
	  /* now compact rmph and rmb to do the remaining integrals */
	  
	  for( i=0, ip=0 ; i<4 ; i++ )
	    {
	      if( ! vremain[i] ) continue ;
	      mvec[ip] = rmb[i] ;
	      for( j=0, jp=0 ; j<4 ; j++ )
		{
		  if( ! vremain[j] ) continue ;
		  mmat[ip][jp] = rmph[i][j] ;
		  jp++ ;
		}
	      ip++ ;
	    }
	  /*
	    To integrate the quadratic form
	    we use the result for M Hermitian
	    Int exp(-(r-ro)M(r-ro)) drn = (pi)^n/2 |M|^-1/2
	    Then  Int exp-(rMr + Br) drn =
	    (pi)^n/2 |M|^-1/2 exp(B(M-1adj)B/4)
	    where B = -2ro M
	  */
	  
	  inv23( ip, mmat, minv ) ;
	  det = det23( ip, mmat ) ;
	  
	  /* adjoint( ip, minv, mmat )   dont need the adjoint */
	  
	  matvec( ip, minv, mvec, xvec ) ;
	  a0 += (dotpro(mvec,xvec))/4. ;
	  countrate *= pow(pi, (double)ip/2.) ;
	  
	  if( a0 > 26. ) a0 = 26. ;
	  if( a0 < -26. ) a0 = -26. ;
	  relint = exp(a0) ;
	  countrate *= relint/sqrt(det) ;
	  Fsqeff = 0 ;
	  *v = -countrate ;
	  retval = 1 ;
	  goto ctrRet ;
	}
    }
  
  /* to get here, no linear condition and icstype not phonon */
  
  countrate *= vol ;
  *v = -countrate ;
  Fsqeff = 0. ;
  /* last case is 1 barn/meV uniform cross-section */
  retval = 1 ;
  
 ctrRet:
  if( gequp ) {
    for( i=0 ; i<3 ; i++ )
      if( gequp[i].coef ) free(gequp[i].coef) ;
    free(gequp) ;
    gequp = NULL ;
  }
  if( rmwl ) getmatrix(rmwl, 3, 0, 0) ;
  if( rmph ) getmatrix(rmph, 4, 0, 0) ;
  if( dmat ) getmatrix(dmat, 3, 0, 0) ;
  if( dinv ) getmatrix(dinv, 3, 0, 0) ;
  if( mmat ) getmatrix(mmat, 3, 0, 0) ;
  if( minv ) getmatrix(minv, 3, 0, 0) ;
  if( Sm ) getmatrix(Sm, 2, 0, 0) ;
  
  rmwl = NULL ;
  rmph = NULL ;
  dmat = NULL ;
  dinv = NULL ;
  mmat = NULL ;
  minv = NULL ;
  Sm = NULL ;
  return retval ;
}


/*
  calc_resm_aux gets rm for final return of proj and intersect and
  stuff for the cross-section plotting and check cs directions
  and gets ki kf in sample coordinates for plotting,
  all at the call point = scan.
  Also calcs global kimag kfmag kisq kfsq qsqr
  N.B. plotting coordinate sys is xuni yuni and using rm
  N.B. qmag and quni globals are calculated at scan
  May 2000 eigen calc moved to resm
*/
int calc_resm_aux(double scan[])
{
  int i, j, k ;
  double omr, alpha, beta ;
  double den, maga, magb ;
  double qpub[3] ;
  double angdum[4], coskfang, sinkfang ;
  /*static double fact = 0.83255461 ; sqrt(ln2) for hwhm */
  double twopi ;
  twopi = TWOPI ;

  for( i=0 ; i<2 ; i++ )
    {
      kixy[i] = 0. ;
      kfxy[i] = 0. ;
      q0[i]   = 0. ;
      cspoint[i] = 0. ;
      cs1diru[i] = 0. ;
      csdire1[i] = 0. ;
      csdire2[i] = 0. ;
      ref1xy[i]  = 0. ;
      ref2xy[i]  = 0. ;
      ref3xy[i]  = 0. ;
      for( j=0 ; j<2 ; j++ ) { sedir[i][j] = 0. ; pedir[i][j] = 0. ; }
    }
  for( i=0 ; i<3 ; i++ )
    {
      axproj[i] = 0. ;
      for( j=0 ; j<2 ; j++ )
	for( k=0 ; k<2 ; k++ )
	  { pr[i][j][k] = 0. ; sc[i][j][k] = 0. ; }
    }
  lor = 0. ;
  cor = 0. ;
  qmag = 0. ;
  stepmag = 0. ;

  /*
    first get the resolution matrix rm 
    at the scan point in the requested coordinate system
    abflag<0 antiCN =0 CN  or >0 ahkl apun  all rhcs
  */
  if( ! calc_resm( scan, iax ) ) return 0 ;
  hklmag = recvec(hkle,quni) ;
  qmag = recvec(scan,quni) ;
  ihklmag = qmag ;
  qsqr = qmag*qmag ;
  stepmag = recvec(step,stepuni) ;


  /* get hkl perp to scan. complete rhc for quni with zuni up */
  crosspro(zuni,quni,qpun) ;
  for( i=0 ; i<3 ; i++ ) qpub[i] = qpun[i] ;
  if( spin[1] == -1 ) negvec(qpub) ;
  /* get hkl for q-perp unit at lower angle omega only for return purposes */

  hklqpun[0] = dotpro(avec,qpub)/twopi ;
  hklqpun[1] = dotpro(bvec,qpub)/twopi ;
  hklqpun[2] = dotpro(cvec,qpub)/twopi ;
  qpmag = 1. ;


  /* aom = scan[3]/DNEUT */
  if( ! kif(ifix,efixd,scan[3],&kimag,&kfmag) ) { scan[3] = 0. ; }
  specangles(scan, angdum) ;
  specanglesd(scan, angls) ;

  /* ttr      = angdum[2] */                        /* in radians */
  omr      = angdum[1] ;                         /* in radians */

  specanglesd(point, angdum) ;
  tth0 = angdum[2] ;

  /*
    corrections for spin[1] removed from angle calc since xtal phkl
    is always lower omega angle compared to ahkl RWE June 1998
    phkl is hkl for apub = unit vector perp to ahkl and towards bhkl
  */
  /*
    hkla already calc unit in set-recip
    for( i=0 ; i<3 ; i++ ) {
    uhkl[i] = ahkl[i]/amag ;
    hkla[i] = uhkl[i] ;
    }
  */

  alpha = -sin(omr) ;
  beta  = cos(omr) ;
  for(i=0 ; i<3 ; i++) hkleI[i] = kimag*(alpha*hkla[i] + beta*phkl[i]);
  hkleI[3] = DNEUT*kimag*kimag ;

  coskfang = (kimag*kimag + kfmag*kfmag - qsqr)/(2.*kimag*kfmag) ;
  sinkfang = sin(acos(coskfang)) ;

  for(i=0 ; i<3 ; i++)
    hkleF[i] = kfmag*((alpha*coskfang + beta*sinkfang)*hkla[i] + 
		      (beta*coskfang - alpha*sinkfang)*phkl[i]) ;
  hkleF[3] = DNEUT*kfmag*kfmag ;


  /*
    now get projections and intersections
    sc pr are for rm proj and intersect with xye, while sedir pedir are
    intersect and project onto Qdir E plane
  */
  /* qdir[0] = 1. ; qdir[1] = 0. ; */
  if( ! proji( rm, axproj, sc, pr ) )
    /*! projidire( rm, qdir, sedir, pedir ) ) */
    {
      setErrmsg(errstr[9]) ;
      return (0) ;
    }

  fwhmeinc = 2.*axproj[2] ;
  fwhmz = 2.*vres ;
  den = sqrt(sc[2][0][0]*sc[2][0][0] + sc[2][1][0]*sc[2][1][0]) ;
  if( den > 0. )
    fwhmb[1] = 2.*fabs(sc[2][1][0]*sc[2][0][1]+sc[2][0][0]*sc[2][1][1])/den ;
  den = sqrt(sc[2][0][1]*sc[2][0][1] + sc[2][1][1]*sc[2][1][1]) ;
  if( den > 0. )
    fwhmb[0] = 2.*fabs(sc[2][1][1]*sc[2][0][0]+sc[2][0][1]*sc[2][1][0])/den ;
  den = sqrt(sc[1][0][1]*sc[1][0][1] + sc[1][1][1]*sc[1][1][1]) ;
  if( den > 0. )
    fwhmb[2] = 2.*fabs(sc[1][1][1]*sc[1][0][0]+sc[1][0][1]*sc[1][1][0])/den ;

  maga = 2.*sqrt(pr[2][0][0]*pr[2][0][0] + pr[2][0][1]*pr[2][0][1]) ;
  magb = 2.*sqrt(pr[2][1][0]*pr[2][1][0] + pr[2][1][1]*pr[2][1][1]) ;
  if( maga > magb ) {
    fwhmp[0] = maga ;
    fwhmp[1] = magb ;
    fwhmp[2] = atan2(pr[2][0][1], pr[2][0][0]) ;
  } else {
    fwhmp[0] = magb ;
    fwhmp[1] = maga ;
    fwhmp[2] = atan2(pr[2][1][1], pr[2][1][0]) ;
  }
  fwhmp[2] *= RADTODEG ;

  lor = 0. ;
  cor = 0. ;
  if( icstyp == 0 )
    lorentz( etok(efixd),&lor,&cor ) ;	/* get Lorentz facs */

  return (1) ;
}

int plotstuff()
{
  int i, j, k, done ;
  double csmag ;

  double hklptuni[3], hkld1uni[3], hkld2uni[3], hkld3uni[3] ;
  double hklpoint[3], hkldire1[3], hkldire2[3], hkldire3[3] ;
  double hklref1[3],  hklref2[3],  hklref3[3] ;
  double hklr1uni[3], hklr2uni[3], hklr3uni[3] ;
  double hkleIuni[3], hkleFuni[3] ;
  double hkleImag,    hkleFmag ;

  double csdir[3][4], fac ;

  /*static double fact = 0.83255461 ;  sqrt(ln2) for hwhm */
  /* double twopi */
  double xyeoff0, xyeoff1 ;
  int nrot ;

  static double **sadj = (double**)0 ;
  static double **s = (double**)0 ;

  if( sadj == (double**)0 )
    {
      sadj = getmatrix(sadj, 0, 3, 3) ;
      s = getmatrix(s, 0, 3, 3) ;
    }

  /* twopi = TWOPI */


  /* setup cross section and ref inplane for plotting in xuni yuni coords */
  for( i=0 ; i<3 ; i++ )
    {
      hklpoint[i] = point[i] ;
      hkldire1[i] = dire1[i] ;
      hkldire2[i] = dire2[i] ;
      hkldire3[i] = dire3[i] ;
      hklref1[i]  = ref1[i] ;
      hklref2[i]  = ref2[i] ;
      hklref3[i]  = ref3[i] ;
    }
  pointmag = recvec(hklpoint,hklptuni) ;
  diremag[0] = recvec(hkldire1,hkld1uni) ;
  diremag[1] = recvec(hkldire2,hkld2uni) ;
  diremag[2] = recvec(hkldire3,hkld3uni) ;
  ref1mag  = recvec(hklref1, hklr1uni) ;
  ref2mag  = recvec(hklref2, hklr2uni) ;
  ref3mag  = recvec(hklref3, hklr3uni) ;

  /* x is ahkl for right scattering and y is ahkl for left */

  hkleImag = recvec(hkleI, hkleIuni) ;
  hkleFmag = recvec(hkleF, hkleFuni) ;
  kixy[0] = hkleImag*dotpro(hkleIuni,xuni) ;
  kixy[1] = hkleImag*dotpro(hkleIuni,yuni) ;
  kfxy[0] = hkleFmag*dotpro(hkleFuni,xuni) ;
  kfxy[1] = hkleFmag*dotpro(hkleFuni,yuni) ;

  q0[0] = qmag*dotpro(quni,xuni) ;
  q0[1] = qmag*dotpro(quni,yuni) ;

  /* use part of cs point in scattering plane */
  cspoint[0] = pointmag*dotpro(hklptuni,xuni) ;
  cspoint[1] = pointmag*dotpro(hklptuni,yuni) ;
  cspoint[2] = pointmag*dotpro(hklptuni,zuni) ;

  /* the q-part of dire1 also defines Q-direction for Q-E 2-plots */
  cs1diru[0] = dotpro(hkld1uni,xuni) ;
  cs1diru[1] = dotpro(hkld1uni,yuni) ;

  csdir[0][0] = diremag[0]*cs1diru[0] ;
  csdir[0][1] = diremag[0]*cs1diru[1] ;
  csdir[0][2] = diremag[0]*dotpro(hkld1uni,zuni) ;
  csdir[0][3] = dire1[3] ;
  csdir[1][0] = diremag[1]*dotpro(hkld2uni,xuni) ;
  csdir[1][1] = diremag[1]*dotpro(hkld2uni,yuni) ;
  csdir[1][2] = diremag[1]*dotpro(hkld2uni,zuni) ;
  csdir[1][3] = dire2[3] ;
  csdir[2][0] = diremag[2]*dotpro(hkld3uni,xuni) ;
  csdir[2][1] = diremag[2]*dotpro(hkld3uni,yuni) ;
  csdir[2][2] = diremag[2]*dotpro(hkld3uni,zuni) ;
  csdir[2][3] = dire3[3] ;
  
  csdire1[0] = csdir[0][0] ;
  csdire1[1] = csdir[0][1] ;
  csdire2[0] = csdir[1][0] ;
  csdire2[1] = csdir[1][1] ;
  
  /* following is to get E=0 intersect scatt plane part of CS */

  if( icstyp == 1 ) /* Bragg line, get in scatt plane part */
    {
      if( diremag[0]>0. && fabs(csdir[0][2]/diremag[0]) > 0.001 )
	{
	  fac = cspoint[2]/csdir[0][2] ;
	  for( i=0 ; i<2 ; i++ )
	    {
	      cspoint[i] -= fac*csdir[0][i] ;
	      csdire1[i] = 0. ;
	    }
	}
      cspoint[3] = 0. ;
      csdire1[3] = 0. ;
    }
  
  if( icstyp == 2 ) /* Bragg plane */
    {
      done = 0 ;
      for( i=0 ; i<2 && done<1 ; i++ )
	{
	  if( diremag[i]>0. && fabs(csdir[i][2]/diremag[i]) > 0.001 )
	    {
	      fac = cspoint[2]/csdir[i][2] ;
	      for( j=0 ; j<2 ; j++ ) cspoint[j] -= fac*csdir[i][j] ;
	      fac = csdir[1-i][2]/csdir[i][2] ;
	      for( j=0 ; j<2 ; j++ )
		{
		  csdire1[j] = csdir[1-i][j] - fac*csdir[i][j] ;
		  csdire2[j] = 0. ;		
		}
	      done = 1 ;
	    }
	}
      cspoint[3] = 0. ;
      csdire1[3] = 0. ;
      csdire2[3] = 0. ;
    }
  
  if( icstyp == 3 )
    {
      csdire1[3] = 0. ;
      csdire2[3] = 0. ;
    }
  
  if( icstyp == 4 )
    {
      cspoint[3] = point[3] ;
      done = 0 ;
      i = 0 ;
      while( i<3 && ! done )
	{
	  if( diremag[i]>0. && fabs(csdir[i][2]/diremag[i]) > 0.001 )
	    {
	      fac = cspoint[2]/csdir[i][2] ;
	      for( j=0 ; j<4 ; j++ ) cspoint[j] -= fac*csdir[i][j] ;
	      j = (i+1) % 3 ;
	      fac = csdir[j][2]/csdir[i][2] ;
	      for( k=0 ; k<4 ; k++ ) 
		csdire1[k] = csdir[j][k] - fac*csdir[i][k] ;
	      j = (i+2) % 3 ;
	      fac = csdir[j][2]/csdir[i][2] ;
	      for( k=0 ; k<4 ; k++ ) 
		csdire2[k] = csdir[j][k] - fac*csdir[i][k] ;
	      done = 1 ;
	    }
	  i++ ;
	}
    }
  
  ref1xy[0]  = ref1mag*dotpro(hklr1uni,xuni) ;
  ref1xy[1]  = ref1mag*dotpro(hklr1uni,yuni) ;
  ref2xy[0]  = ref2mag*dotpro(hklr2uni,xuni) ;
  ref2xy[1]  = ref2mag*dotpro(hklr2uni,yuni) ;
  ref3xy[0]  = ref3mag*dotpro(hklr3uni,xuni) ;
  ref3xy[1]  = ref3mag*dotpro(hklr3uni,yuni) ;
  
  if(reduced) {
    kixy[0] /= xmag ;
    kixy[1] /= ymag ;
    kfxy[0] /= xmag ;
    kfxy[1] /= ymag ;
    q0[0] /= xmag ;
    q0[1] /= ymag ;
    cspoint[0] /= xmag ;
    cspoint[1] /= ymag ;
    cs1diru[0] /= xmag ;
    cs1diru[1] /= ymag ;
    csmag = sqrt(cs1diru[0]*cs1diru[0] + cs1diru[1]*cs1diru[1]) ;
    cs1diru[0] /= csmag ;
    cs1diru[1] /= csmag ;
    
    csdire1[0] /= xmag ;
    csdire1[1] /= ymag ;
    csdire2[0] /= xmag ;
    csdire2[1] /= ymag ;
    
    ref1xy[0]  /= xmag ;
    ref1xy[1]  /= ymag ;
    ref2xy[0]  /= xmag ;
    ref2xy[1]  /= ymag ;
    ref3xy[0]  /= xmag ;
    ref3xy[1]  /= ymag ;
  }

    
  /* now transform rmcn to x and y coordinates for plotting */
  s[0][0] = -dotpro(quni,xuni) ;
  s[0][1] = -dotpro(quni,yuni) ;
  s[0][2] = 0. ;
  s[1][0] = -dotpro(qpun,xuni) ;
  s[1][1] = -dotpro(qpun,yuni) ;
  s[1][2] = 0. ;
  s[2][0] = 0. ;
  s[2][1] = 0. ;
  s[2][2] = 1. ;
  for( i=0 ; i<3 ; i++ )
    for( j=0 ; j<3 ; j++ )
      { sadj[i][j] = s[i][j] ; }
  sadj[0][1] = s[1][0] ;
  sadj[1][0] = s[0][1] ;
  matmat3(sadj, rmcn, rmp) ;
  matmat3(rmp, s, rmplot) ;


  xyeoff0 = xyeoff[0] ;
  xyeoff1 = xyeoff[1] ;
  xyeoffp[0] = s[0][0]*xyeoff0 + s[1][0]*xyeoff1 ;
  xyeoffp[1] = s[0][1]*xyeoff0 + s[1][1]*xyeoff1 ;
  
  if(reduced) {
    s[0][0] = xmag ; s[0][1] = 0. ;
    s[1][0] = 0. ;   s[1][1] = ymag ;
    matmat3(s, rmplot, rmp) ;
    matmat3(rmp, s, rmplot) ;
    xyeoffp[0] /= xmag ;
    xyeoffp[1] /= ymag ;
  }
  q0[0] += xyeoffp[0] ;
  q0[1] += xyeoffp[1] ;

  if( ! proji( rmplot, axprojp, scp, prp ) ||
      ! projidire( rmplot, q0, sedirp, pedirp ) )
    {
      setErrmsg(errstr[9]) ;
      getmatrix(sadj, 3, 0, 0) ;
      getmatrix(s, 3, 0, 0) ;
      sadj = NULL ;
      s = NULL ;
      return (0) ;
    }

  /*
    we need the eigenvalues and norm-eigenvectors to plot
    eigenvecs is the rotation matrix for GL
  */

  jacobi( rmplot, 3, eigenvals, eigenvecs, &nrot ) ;
  getmatrix(sadj, 3, 0, 0) ;
  getmatrix(s, 3, 0, 0) ;
  sadj = NULL ;
  s = NULL ;
  return(1) ;
}

int calc_resm( double scan[], int abflag )
{
  /*
    calc_resm
    calculates the res matrix, rm,  in C&N or antiC&N
    or if abflag calcs rmw in ahkl ahklperp as used by width calc
  */
  int i, j ;
  double qm, eI, eF ;
  static double **sadj = (double**)0 ;
  static double **s = (double**)0 ;



  //static double fact = 0.83255461 ;  /* sqrt(ln2) for hwhm */

  if( sadj == (double**)0 )
    {
      sadj = getmatrix(sadj, 0, 3, 3) ;
      s = getmatrix(s, 0, 3, 3) ;
    }
  qm = recvec(scan,quni) ;
  if( ifixmon ) { eI = efixd ; eF = efixd - *(scan+3) ; }
  else { eF = efixd ; eI = efixd + *(scan+3) ; }

  if( ! resm( eI, eF, qm ) ) {
    getmatrix(sadj, 3, 0, 0) ;
    getmatrix(s, 3, 0, 0) ;
    sadj = NULL ;
    s = NULL ;
    return 0 ;
  }
  for( i=0 ; i<3 ; i++ )
    for( j=0 ; j<3 ; j++ ) rmcn[i][j] = rm[i][j] ;
 
  /*
    resm returns rm in C&N coordinates with Qx antiparallel to Q and
    Qy completing rhc with z up.
    abflag is used to flag coordinate system change
    to ahkl and bhkl-perp (A-1 units) completing rhc with z-up.
    Recall that bhkl is to lower angle in omega compared to ahkl, so
    for left scattering (-1) bhkl-perp is x-axis, ahkl is y-axis
    for rght scattering ( 1) ahkl      is x-axis, bhkl-perp is y-axis
    */

  if( abflag < 0 ) {
    rm[0][2] = -rm[0][2] ; rm[2][0] = -rm[2][0] ;
    rm[1][2] = -rm[1][2] ; rm[2][1] = -rm[2][1] ;
  }

  if( abflag <= 0 ) {
    getmatrix(sadj, 3, 0, 0) ;
    getmatrix(s, 3, 0, 0) ;
    sadj = NULL ;
    s = NULL ;
    return (1) ;
  }
  /*
    transform to ahkl <-> bhkl-perp (apub) directions for plotting.
    and to x=ahkl y=apun for width calculation in rhc

    radj RM r = s   where r = (Q, Qperp, E)
    to transform to x y the full transformation matrix S is
    Qhat . xhat      Qhat . yhat        0
    Qperphat . xhat  Qperphat . yhat    0
    0                0                      1
    so (Q, Qperp) = S x (x, y) regardless if x, y is rhcsys
    or Sadj x (Q, Qperp) = (x, y)
    so radj (R = S Sadj R S Sadj) r = s = (Sadj r)adj (Sadj R S) (Sadj r)
    this is the equ in x,y coordinates
    If we want to go to reduced units so that
    (xred, yred) = F Sadj x (Q, Qperp) 
    with F11 = 1/xmag and F22 = 1/ymag
    the new ellipsoid equ is F-1 (Sadj R S) F-1
  */
  

  /*
    get qpun, perp to quni completing rhc with x and y axes negative of C&N
  */
  crosspro(zuni,quni,qpun) ;

  /*
    first transform rm to auni and apun 
    rh-coordinates for width calc with rmw
  */
  s[0][0] = -dotpro(quni,auni) ;
  s[0][1] = -dotpro(quni,apun) ;
  s[0][2] = 0. ;
  s[1][0] = -dotpro(qpun,auni) ;
  s[1][1] = -dotpro(qpun,apun) ;
  s[1][2] = 0. ;
  s[2][0] = 0. ;
  s[2][1] = 0. ;
  s[2][2] = 1. ;
  for( i=0 ; i<3 ; i++ )
    for( j=0 ; j<3 ; j++ )
      { sadj[i][j] = s[i][j] ; }
  sadj[0][1] = s[1][0] ;
  sadj[1][0] = s[0][1] ;
  matmat3(sadj, rm, rmp) ;
  matmat3(rmp, s, rmw) ;

  /* also transform xyeoff for width calc */
  xyeoffw[0] = s[0][0]*xyeoff[0] + s[1][0]*xyeoff[1] ;
  xyeoffw[1] = s[0][1]*xyeoff[0] + s[1][1]*xyeoff[1] ;
  xyeoffw[2] = xyeoff[2] ;

  /*
    get hkl for xyeoffw
    auni = hkla[0]ast + hkla[1]bst + hkla[2]cst
    apun = hklp[0]ast + hklp[1]bst + hklp[2]cst
    zuni = hklz[0]ast + hklz[1]bst + hklz[2]cst
    xyeoffw = xyeoffw[0]*auni + xyeoffw[1]*apun + qzoff*zuni
            = (xyeoffw[0]*hkla[0] + xyeoffw[1]*hklp[0] + qzoff*hklz[0])ast
            + (xyeoffw[0]*hkla[1] + xyeoffw[1]*hklp[1] + qzoff*hklz[1])bst
	    + (xyeoffw[0]*hkla[2] + xyeoffw[1]*hklp[2] + qzoff*hklz[2])cst
  */
  for( i=0 ; i<3 ; i++ ) {
    hkloff[i] = xyeoffw[0]*hkla[i] + xyeoffw[1]*hklp[i] + qzoff*hklz[i] ;
    hkleff[i] = scan[i] + hkloff[i] ;
  }
  hkleff[3] = scan[3] + xyeoffw[2] ;

  /* copy rmw to rm and xyeoffw to xyeoff  as requested by abflag > 0 */

  for( i=0 ; i<3 ; i++ ) {
    xyeoff[i] = xyeoffw[i] ;
    for( j=0 ; j<3 ; j++ ) rm[i][j] = rmw[i][j] ;
  }
  getmatrix(sadj, 3, 0, 0) ;
  getmatrix(s, 3, 0, 0) ;
  sadj = NULL ;
  s = NULL ;
  return (1) ;
}

static int calc_resmX( double scan[], int abflag )
{
  /*
    calc_resmX
    calculates the res matrix, rm,  in C&N or antiC&N
    or if abflag calcs rmw in ahkl ahklperp as used by width calc
    This version uses global Ki and Kf instead of ifixmon and efixd
  */
  int i, j ;
  double qm, eI, eF ;
  static double **sadj = (double**)0 ;
  static double **s = (double**)0 ;
  static double Dn = 2.072141789 ;

  //static double fact = 0.83255461 ;  /* sqrt(ln2) for hwhm */

  if( sadj == (double**)0 )
    {
      sadj = getmatrix(sadj, 0, 3, 3) ;
      s = getmatrix(s, 0, 3, 3) ;
    }
  qm = recvec(scan,quni) ;

  eI = Dn*Ki*Ki ;
  eF = Dn*Kf*Kf ;
  /* NB we dont check consistency of eI - eF and scan[3] */

  if( ! resm( eI, eF, qm ) ) {
    getmatrix(sadj, 3, 0, 0) ;
    getmatrix(s, 3, 0, 0) ;
    sadj = NULL ;
    s = NULL ;
    return 0 ;
  }
  for( i=0 ; i<3 ; i++ )
    for( j=0 ; j<3 ; j++ ) rmcn[i][j] = rm[i][j] ;
 
  /*
    resm returns rm in C&N coordinates with Qx antiparallel to Q and
    Qy completing rhc with z up.
    abflag is used to flag coordinate system change
    to ahkl and bhkl-perp (A-1 units) completing rhc with z-up.
    Recall that bhkl is to lower angle in omega compared to ahkl, so
    for left scattering (-1) bhkl-perp is x-axis, ahkl is y-axis
    for rght scattering ( 1) ahkl      is x-axis, bhkl-perp is y-axis
    */

  if( abflag < 0 ) {
    rm[0][2] = -rm[0][2] ; rm[2][0] = -rm[2][0] ;
    rm[1][2] = -rm[1][2] ; rm[2][1] = -rm[2][1] ;
  }

  if( abflag <= 0 ) {
    getmatrix(sadj, 3, 0, 0) ;
    getmatrix(s, 3, 0, 0) ;
    sadj = NULL ;
    s = NULL ;
    return (1) ;
  }
  /*
    transform to ahkl <-> bhkl-perp (apub) directions for plotting.
    and to x=ahkl y=apun for width calculation in rhc

    radj RM r = s   where r = (Q, Qperp, E)
    to transform to x y the full transformation matrix S is
    Qhat . xhat      Qhat . yhat        0
    Qperphat . xhat  Qperphat . yhat    0
    0                0                      1
    so (Q, Qperp) = S x (x, y) regardless if x, y is rhcsys
    or Sadj x (Q, Qperp) = (x, y)
    so radj (R = S Sadj R S Sadj) r = s = (Sadj r)adj (Sadj R S) (Sadj r)
    this is the equ in x,y coordinates
    If we want to go to reduced units so that
    (xred, yred) = F Sadj x (Q, Qperp) 
    with F11 = 1/xmag and F22 = 1/ymag
    the new ellipsoid equ is F-1 (Sadj R S) F-1
  */
  

  /*
    get qpun, perp to quni completing rhc with x and y axes negative of C&N
  */
  crosspro(zuni,quni,qpun) ;

  /*
    first transform rm to auni and apun 
    rh-coordinates for width calc with rmw
  */
  s[0][0] = -dotpro(quni,auni) ;
  s[0][1] = -dotpro(quni,apun) ;
  s[0][2] = 0. ;
  s[1][0] = -dotpro(qpun,auni) ;
  s[1][1] = -dotpro(qpun,apun) ;
  s[1][2] = 0. ;
  s[2][0] = 0. ;
  s[2][1] = 0. ;
  s[2][2] = 1. ;
  for( i=0 ; i<3 ; i++ )
    for( j=0 ; j<3 ; j++ )
      { sadj[i][j] = s[i][j] ; }
  sadj[0][1] = s[1][0] ;
  sadj[1][0] = s[0][1] ;
  matmat3(sadj, rm, rmp) ;
  matmat3(rmp, s, rmw) ;

  /* also transform xyeoff for width calc */
  xyeoffw[0] = s[0][0]*xyeoff[0] + s[1][0]*xyeoff[1] ;
  xyeoffw[1] = s[0][1]*xyeoff[0] + s[1][1]*xyeoff[1] ;
  xyeoffw[2] = xyeoff[2] ;

  /*
    get hkl for xyeoffw
    auni = hkla[0]ast + hkla[1]bst + hkla[2]cst
    apun = hklp[0]ast + hklp[1]bst + hklp[2]cst
    zuni = hklz[0]ast + hklz[1]bst + hklz[2]cst
    xyeoffw = xyeoffw[0]*auni + xyeoffw[1]*apun + qzoff*zuni
            = (xyeoffw[0]*hkla[0] + xyeoffw[1]*hklp[0] + qzoff*hklz[0])ast
            + (xyeoffw[0]*hkla[1] + xyeoffw[1]*hklp[1] + qzoff*hklz[1])bst
	    + (xyeoffw[0]*hkla[2] + xyeoffw[1]*hklp[2] + qzoff*hklz[2])cst
  */
  for( i=0 ; i<3 ; i++ ) {
    hkloff[i] = xyeoffw[0]*hkla[i] + xyeoffw[1]*hklp[i] + qzoff*hklz[i] ;
    hkleff[i] = scan[i] + hkloff[i] ;
  }
  hkleff[3] = scan[3] + xyeoffw[2] ;

  /* copy rmw to rm and xyeoffw to xyeoff  as requested by abflag > 0 */

  for( i=0 ; i<3 ; i++ ) {
    xyeoff[i] = xyeoffw[i] ;
    for( j=0 ; j<3 ; j++ ) rm[i][j] = rmw[i][j] ;
  }
  getmatrix(sadj, 3, 0, 0) ;
  getmatrix(s, 3, 0, 0) ;
  sadj = NULL ;
  s = NULL ;
  return (1) ;
}


/*
  note that rmw matrix is used in rmvec for width calcs
  Sep 99 RWE change vectors so 0,1,2 are Qx Qy Qz and 3=E
*/
void rmvec( double vecin[], double vecout[] )
{
	vecout[0] = 
	  rmw[0][0]*vecin[0] + rmw[0][1]*vecin[1] + rmw[0][2]*vecin[3] ;
	vecout[1] = 
	  rmw[1][0]*vecin[0] + rmw[1][1]*vecin[1] + rmw[1][2]*vecin[3] ;
	vecout[2] = rmv*vecin[2] ;
	vecout[3] = 
	  rmw[2][0]*vecin[0] + rmw[2][1]*vecin[1] + rmw[2][2]*vecin[3] ;
}

double vrv( double vec1[], double vec2[] )
{
	double v[4] ;
	rmvec( vec2, v ) ;
	return (dotpro4( v, vec1 )) ;
}

void rmqvec( double vecin[], double vecout[] )
{
	int i, j ;
	for( i=0 ; i<2 ; i++ )
		for( j=0, vecout[i] = 0. ; j<2 ; j++ )
			vecout[i] += rmw[i][j]*vecin[j] ;
}

double vrvq( double vec1[], double vec2[] )
{
	double v[2] ;
	rmqvec( vec2, v ) ;
	return (v[0]*vec1[0] + v[1]*vec1[1] + rmv*vec1[2]*vec2[2]) ;
}






#define PRINTERR printf("%s\n",errmsg)
#define ERR { PRINTERR ; continue ; }


int set_fix()
{
  if	( ifixmon == 1 )	ifix = 2 ;
  else	ifix = 3 ;
  return (1) ;
}

int lowercase( char *bstr )
{
  char *cp ;
  cp = bstr ;
  while( *cp != '\0' ) { *cp = tolower(*cp) ; cp++ ; }
  return (1) ;
}

/* remember that parser converts all data to lowercase */
int set_cs()
{
  lowercase(cstyp) ;
  if ( strstr(cstyp, "point") != NULL ) icstyp = 0 ;         /* Bragg point */
  else if ( strstr(cstyp, "plane") != NULL ) icstyp = 2 ;
  else if ( strstr(cstyp, "line")  != NULL && strstr(cstyp, "disp") == NULL )
    icstyp = 1 ;
  else if ( strstr(cstyp, "q") != NULL && strstr(cstyp, "fix") != NULL )
    icstyp = 3 ;  /* fix e q-ind */
  else if ( strstr(cstyp, "linear") != NULL || strstr(cstyp, "phon") != NULL )
    icstyp = 4 ;   /* linear dispersion */
  else icstyp = 5 ;            /* energy and Q independent */
  return (1) ;
}

int set_xaxis()
{
  iax = 0 ;
  if ( *xaxistyp == 'a' ) iax = 1 ;
  else if( *xaxistyp == '-' ) iax = -1 ;
  return (1) ;
}

int set_spin()
{
  if( sscanf(aspin, "%d %d %d", spin, spin+1, spin+2) == 3 ) return (1) ;
  spin[0] = 1 ; if(*aspin=='r' || *aspin== 'R') spin[0] = -1 ;
  spin[1] = -1 ; if(*(aspin+1)=='l' || *(aspin+1)== 'L') spin[1] = 1 ;
  spin[2] = 1 ; if(*(aspin+2)=='r' || *(aspin+2)== 'R') spin[2] = -1 ;
  scattsideright = 0 ;
  if( spin[1] == -1 ) scattsideright = 1 ;
  return (1) ;
}


char *strcpyend( char *start, char *string )
{
  int nc ;
  if( start == (char*)0 || string == (char*)0 ) return (char*)0 ;
  nc = strlen(string) ;
  strcpy(start, string) ;
  return (start + nc) ;
}

static double xofscan( double scan[] )
{
  int i, nx ;
  double xsum ;
  /*
    for each nonzero step coordinate find x and average
  */
  nx = 0 ;
  xsum = 0. ;
  for( i=0 ; i<4 ; i++ ) {
    if( fabs(step[i]) > 0. ) {
      xsum += (scan[i] - hkle[i])/step[i] ;
      nx++ ;
    }
  }
  if( nx < 1 ) return 0. ;
  return (xsum/nx) ;
}
void Res3_SetDopt(double h, double k, double l, double e)
{
  dopt[0] = h ;
  dopt[1] = k ;
  dopt[2] = l ;
  dopt[3] = e ;
}

static void dores3calc()
{
  set_fix() ;
  set_spin() ;
  set_cs() ;
  set_xaxis() ;
  set_recip() ;
  //printf("dores3calc: about to call width_calc\n") ;
  width_calc() ;
  efixd = efix + xofscan(dopt)*defix ;
  hklecountrate = countrates(hklei) ;
  efixd = efix + xofscan(peak)*defix ;
  maxcountrate = countrates(peak) ;
  /* countrates also does extinction so do this at peak */
 
  /* do all the plotting stuff at dopt */
  efixd = efix + xofscan(dopt)*defix ;

  calc_resm_aux( dopt ) ;
  plotstuff() ;
  /* calc standard C&N matrix and save rm in rmcn */
  /* calc_resm( dopt, 0 ) ; */
  /* skip aux stuff and put in separate command */
}      
static void dores3rm()
{
  int i, j ;

  set_fix() ;
  set_spin() ;
  set_cs() ;
  set_xaxis() ;
  set_recip() ;
  /* width_calc() ; */
  /* hklecountrate = countrates(dopt) ; */
  /* maxcountrate = countrates(peak) ; */
  /* countrates also does extinction so do this at peak */
  /* calc standard C&N matrix and save rm in rmcn */
  calc_resm( dopt, 0 ) ;
  for( i=0 ; i<3 ; i++ )
    for( j=0 ; j<3 ; j++ ) rmcn[i][j] = rm[i][j] ;
  
  /* do all the plotting stuff at dopt */
  /* calc_resm_aux( dopt ) ; */
}      
static int dores3rmab()
{
  /* for use by multixtal */
  return calc_resmX( dopt, 1 ) ;
}

double ***makeN22p(double M3[][2][2], int n)
{
  /* make a 3pointer that can be indexed as the static array M3 */
  int i, j ;
  double ***m3 ;
  m3 = (double***)malloc(n*sizeof(double**)) ;
  for( i=0 ; i<n ; i++ ) {
    m3[i] = (double**)malloc(2*sizeof(double*)) ;
    for( j=0 ; j<2 ; j++ ) {
      m3[i][j] = (double*)&M3[i][j] ;
    }
  }  
  return m3 ;
}
void freeN22p(double ***m3, int n)
{
  int i ;
  if( m3 == NULL ) return ;
  for( i=0 ; i<n ; i++ ) if(m3[i] != NULL) free(m3[i]) ;
  free(m3) ;
}
double **makeN3p(double M2[][3], int n)
{
  /* make a 2pointer than can be indexed as the static array M2 */
  int i ;
  double **m2 ;

  m2 = (double**)malloc(n*sizeof(double*)) ;
  for( i=0 ; i<n ; i++ ) {
    m2[i] = (double*)&M2[i] ;
  }
  return m2 ;
}
void freeN3p(double **m2)
{
  if( m2 ) free(m2) ;
}

double **makeN2p(double M2[][2], int n)
{
  /* make a 2pointer than can be indexed as the static array M2 */
  int i ;
  double **m2 ;

  m2 = (double**)malloc(n*sizeof(double*)) ;
  for( i=0 ; i<n ; i++ ) {
    m2[i] = (double*)&M2[i] ;
  }
  return m2 ;
}
void freeN2p(double **m2)
{
  if( m2 ) free(m2) ;
}

/*
  startup inits dictionary to store input data.
*/


/*
  We are using static matrices so we can use C compiler which requires
  constant addresses when initializing
  The make2p ... functions above are for creating an equivalent
  pointer which can be indexed like the static array
  but can also be used in function calls with pointers.
*/


/*
typedef struct {
  double latt[3] ;
  double angl[3] ;
  double ahkl[3] ;
  double bhkl[3] ;
  double hcol[4] ;
  double hgeo[4] ;
  double vcol[4] ;
  double vgeo[4] ;
  double mosa[3] ;
  double vmos[3] ;
  double dsps[2] ;
  double distance[4] ;
  double efix, defix, lamb ;
  int reactormode, samplemode, detectormode ;
  int offseton ;
  int ifixmon ;
  char aspin[4] ;
  char xaxistyp[8] ;
  char *monoprog ;
  char *analprog ;
  char *datasrc ;
  double rxtotalflux ;
  double epeakflux ;
  double monArea ;
  double moneff ;
  double mref ;
  double aref ;
} Conv3_resinfo ;
*/

void Res3_loadResinfo(Conv3_resinfo *respt)
{
  int i ;
  for( i=0 ; i<3 ; i++ ) respt->latt[i] = latt[i] ;
  for( i=0 ; i<3 ; i++ ) respt->angl[i] = angl[i] ;
  for( i=0 ; i<3 ; i++ ) respt->ahkl[i] = ahkl[i] ;
  for( i=0 ; i<3 ; i++ ) respt->bhkl[i] = bhkl[i] ;
  for( i=0 ; i<4 ; i++ ) respt->hcol[i] = hcol[i] ;
  for( i=0 ; i<4 ; i++ ) respt->hgeo[i] = geoalpha[i] ;
  for( i=0 ; i<4 ; i++ ) respt->vcol[i] = vcol[i] ;
  for( i=0 ; i<4 ; i++ ) respt->vgeo[i] = geobeta[i] ;
  for( i=0 ; i<3 ; i++ ) respt->mos[i] = mosa[i] ;
  for( i=0 ; i<3 ; i++ ) respt->vmos[i] = vmos[i] ;
  for( i=0 ; i<2 ; i++ ) respt->dsps[i] = dsps[i] ;
  for( i=0 ; i<4 ; i++ ) respt->distance[i] = distance[i] ;
  respt->efix = efix ;
  respt->defix = defix ;
  respt->lambda = lamb ;
  respt->angleA = angleA ;
  respt->chiA = chiA ;
  respt->chiB = chiB ;
  respt->monomode = monomode ;
  respt->analmode = analmode ;
  respt->geomode = geomode ;
  respt->offseton = offseton ;
  respt->fixm = ifixmon ;
  respt->scattsideright = scattsideright ;
  for( i=0 ; i<32 ; i++ ) respt->aspin[i] = aspin[i] ; 
  for( i=0 ; i<8 ; i++ ) respt->xaxistyp[i] = xaxistyp[i] ;

  if( respt->monoprog != NULL ) {
    free(respt->monoprog) ; respt->monoprog = NULL ;
  }
  if( monoprog == NULL || strlen(monoprog) < 1 ) respt->monoprog = NULL ;
  else respt->monoprog = strdup(monoprog) ;

  if( respt->analprog != NULL ) {
    free(respt->analprog) ; respt->analprog = NULL ;
  }
  if( analprog == NULL || strlen(analprog) < 1 ) respt->analprog = NULL ;
  else respt->analprog = strdup(analprog) ;

  if( respt->datasrc != NULL ) {
    free(respt->datasrc) ; respt->datasrc = NULL ;
  }
  if( datasrc == NULL || strlen(datasrc) < 1 ) respt->datasrc = NULL ;
  else respt->datasrc = strdup(datasrc) ;

  respt->rxtotalflux = rxtotalflux ;
  respt->epeakflux = epeakflux ;
  respt->monArea = monArea ;
  respt->moneff = moneff ;
  respt->mref = mref ;
  respt->aref = aref ;
  for( i=0 ; i<4 ; i++ ) {
    respt->Imax[i] = Imax[i] ;
    respt->Imin[i] = Imin[i] ;
    respt->Irel[i] = Irel[i] ;
    respt->Iabs[i] = Iabs[i] ;
  }
  respt->resprobmin = resprobmin ;
  respt->npts = npts ;
  respt->flds[0] = Temp ;
  respt->flds[1] = Hfld ;
  respt->flds[2] = PolZ ;
  respt->flds[3] = PolX ;
  respt->flds[4] = PolY ;
}
void Res3_unloadResinfo(Conv3_resinfo *respt)
{
  int i ;
  for( i=0 ; i<3 ; i++ ) latt[i] = respt->latt[i] ;
  for( i=0 ; i<3 ; i++ ) angl[i] = respt->angl[i] ;
  for( i=0 ; i<3 ; i++ ) ahkl[i] = respt->ahkl[i] ;
  for( i=0 ; i<3 ; i++ ) bhkl[i] = respt->bhkl[i] ;
  for( i=0 ; i<4 ; i++ ) hcol[i] = respt->hcol[i] ;
  for( i=0 ; i<4 ; i++ ) geoalpha[i] = respt->hgeo[i] ;
  for( i=0 ; i<4 ; i++ ) vcol[i] = respt->vcol[i] ;
  for( i=0 ; i<4 ; i++ ) geobeta[i] = respt->vgeo[i] ;
  for( i=0 ; i<3 ; i++ ) mosa[i] = respt->mos[i] ;
  for( i=0 ; i<3 ; i++ ) vmos[i] = respt->vmos[i] ;
  for( i=0 ; i<2 ; i++ ) dsps[i] = respt->dsps[i] ;
  for( i=0 ; i<4 ; i++ ) distance[i] = respt->distance[i] ;
  efix = respt->efix ;
  defix = respt->defix ;
  lamb = respt->lambda ;
  angleA = respt->angleA ;
  chiA = respt->chiA ;
  chiB = respt->chiB ;
  monomode = respt->monomode ;
  analmode = respt->analmode ;
  geomode = respt->geomode ;
  offseton = respt->offseton ;
  ifixmon = respt->fixm ;
  scattsideright = respt->scattsideright ;
  for( i=0 ; i<32 ; i++ ) aspin[i] = respt->aspin[i] ; 
  for( i=0 ; i<8 ; i++ ) xaxistyp[i] = respt->xaxistyp[i] ;

  if( monoprog != NULL ) {
    free(monoprog) ; monoprog = NULL ;
  }
  if( respt->monoprog == NULL || strlen(respt->monoprog) < 1 )
    monoprog = NULL ;
  else monoprog = strdup(respt->monoprog) ;

  if( analprog != NULL ) {
    free(analprog) ; analprog = NULL ;
  }
  if( respt->analprog == NULL || strlen(respt->analprog) < 1 )
    analprog = NULL ;
  else analprog = strdup(respt->analprog) ;

  if( datasrc != NULL ) {
    free(datasrc) ; datasrc = NULL ;
  }
  if( respt->datasrc == NULL || strlen(respt->datasrc) < 1 )
    datasrc = NULL ;
  else datasrc = strdup(respt->datasrc) ;

  rxtotalflux = respt->rxtotalflux ;
  epeakflux = respt->epeakflux ;
  monArea = respt->monArea ;
  moneff = respt->moneff ;
  mref = respt->mref ;
  aref = respt->aref ;
  for( i=0 ; i<4 ; i++ ) {
    Imax[i] = respt->Imax[i] ;
    Imin[i] = respt->Imin[i] ;
    Irel[i] = respt->Irel[i] ;
    Iabs[i] = respt->Iabs[i] ;
  }
  resprobmin = respt->resprobmin ;
  npts = respt->npts ;
  Temp = respt->flds[0] ;
  Hfld = respt->flds[1] ;
  PolZ = respt->flds[2] ;
  PolX = respt->flds[3] ;
  PolY = respt->flds[4] ;
}


static void startup()
{


  /*
    dictionary  alt
    LATT#       latt*[#] a [b c angBC angAC angAB]
    ANGL#       angl[#]
    AHKL#       ahkl
    BHKL#       bhkl
                orient* OR xstal*  a b c X angBC angAC angAB
                  sets  LATT1 LATT2 LATT3 X ANGL1 ANGL2 ANGL3
    HCOL#       coll*[#]
                hcol*[#]
    VCOL#       vcol*[#]
    GHCO#       geo*
                hgeo*
    GVCO#       vgeo*
    MOSAM       mosm*
    MOSAS       moss*
    MOSAA       mosa*
                mos mosa moss mosa
    VMOSM       vmosm*
    VMOSS       vmoss*
    VMOSA       vmosa*
                vmos vmosm vmoss vmosa
    DSPM        dspm*
                dm
		mdsp*
    DSPA        dspa*
                da
		adsp*
    DRM         dist*1
    DMS         dist*2
    DSA         dist*3
    DAD         dist*4
    ANGLA       angleA*
    CHIA        chi*
                chiA*
    CHIB        chiB*
    LR          spin*
                aspin*
		rl*
		lr*
		scatt*

    EFIX        efix*
                emfix* also sets FIXM
		eafix* also sets FIXM
    DEFIX       defix*
    LAMB        lamb*

    FIXM        fixm*

                mono* prog  special set alloc
                anal* prog  special set alloc

    FLUX        flux
    epeak*
    mref*
    aref*
    ref*
    irel*
    iabs*
    resprob*
    prob*
    offset*
    reactor*
    sample*
    detect*
  */

  /*
  rmp = getmatrix(rmp, 0, 3, 3) ;
  rmw = getmatrix(rmw, 0, 3, 3) ;
  rm44 = getmatrix(rm44, 0, 3, 3) ;
  rmq = getmatrix(rmq, 0, 3, 3) ;

  rm = getmatrix(rm, 0, 3, 3) ;
  rmcn = getmatrix(rmcn, 0, 3, 3) ;

  pr = getmatrix3(pr, 0, 0, 3, 2, 2) ;
  sc = getmatrix3(sc, 0, 0, 3, 2, 2) ;
  sedir = getmatrix(sedir, 0, 2, 2) ;
  pedir = getmatrix(pedir, 0, 2, 2) ;

  eigenvecs = getmatrix(eigenvecs, 0, 3, 3) ;
  */



  rmp = makeN3p(rmpM, 3) ;
  rmplot = makeN3p(rmplotM, 3) ;
  rmw = makeN3p(rmwM, 3) ;
  /* rm44 = makeN3p(rm44M, 3) */
  rmq = makeN3p(rmqM, 3) ;
  rm = makeN3p(rmM, 3) ;
  rmcn = makeN3p(rmcnM, 3) ;
  sedir = makeN2p(sedirM, 2) ;
  pedir = makeN2p(pedirM, 2) ;
  sedirp = makeN2p(sedirpM, 2) ;
  pedirp = makeN2p(pedirpM, 2) ;
  eigenvecs = makeN3p(eigenvecsM, 3) ;

  pr = makeN22p(prM, 3) ;
  sc = makeN22p(scM, 3) ;
  prp = makeN22p(prpM, 3) ;
  scp = makeN22p(scpM, 3) ;

  /*
    equ[0].coef = (double*)malloc(4*sizeof(double)) ;
    equ[1].coef = (double*)malloc(4*sizeof(double)) ;
  */
  
  globalinput = &inputentries[0] ;
  globaloutput = &outputentries[0] ;

  globalinputdict = &inputdict ;
  globaloutputdict = &outputdict ;

  /* load the input and output cmnds into dictionarys */
  dictstartup( globalinputdict, inputentries ) ;
  dictstartup( globaloutputdict, outputentries ) ;
}

void Res3_shutdown()
{
  freeN3p(rmp) ;
  freeN3p(rmplot) ;
  freeN3p(rmw) ;
  freeN3p(rmq) ;
  freeN3p(rm) ;
  freeN3p(rmcn) ;
  freeN2p(sedir) ;
  freeN2p(pedir) ;
  freeN2p(sedirp) ;
  freeN2p(pedirp) ;
  freeN3p(eigenvecs) ;

  freeN22p(pr, 3) ;
  freeN22p(sc, 3) ;
  freeN22p(prp, 3) ;
  freeN22p(scp, 3) ;

  /* kill dictionarys */
  dictkill( globalinputdict, inputentries ) ;
  dictkill( globaloutputdict, outputentries ) ;
}










