#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _rcosmo_car2sph(SEXP);
extern SEXP _rcosmo_covCMB_internal1(SEXP, SEXP);
extern SEXP _rcosmo_covCMB_internal2(SEXP, SEXP);
extern SEXP _rcosmo_maxDist_internal(SEXP);
extern SEXP _rcosmo_minDist_internal(SEXP, SEXP);
extern SEXP _rcosmo_mkpix2xyC(SEXP);
extern SEXP _rcosmo_nest2ring(SEXP, SEXP);
extern SEXP _rcosmo_pix2coords_internal(SEXP, SEXP, SEXP, SEXP);
extern SEXP _rcosmo_pointInConvexPolygon(SEXP, SEXP);
extern SEXP _rcosmo_pointInConvexPolygonHP(SEXP, SEXP, SEXP, SEXP);
extern SEXP _rcosmo_pointInDisc(SEXP, SEXP);
extern SEXP _rcosmo_pointInDiscHP(SEXP, SEXP, SEXP, SEXP);
extern SEXP _rcosmo_sph2car(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_rcosmo_car2sph",                (DL_FUNC) &_rcosmo_car2sph,                1},
  {"_rcosmo_covCMB_internal1",       (DL_FUNC) &_rcosmo_covCMB_internal1,       2},
  {"_rcosmo_covCMB_internal2",       (DL_FUNC) &_rcosmo_covCMB_internal2,       2},
  {"_rcosmo_maxDist_internal",       (DL_FUNC) &_rcosmo_maxDist_internal,       1},
  {"_rcosmo_minDist_internal",       (DL_FUNC) &_rcosmo_minDist_internal,       2},
  {"_rcosmo_mkpix2xyC",              (DL_FUNC) &_rcosmo_mkpix2xyC,              1},
  {"_rcosmo_nest2ring",              (DL_FUNC) &_rcosmo_nest2ring,              2},
  {"_rcosmo_pix2coords_internal",    (DL_FUNC) &_rcosmo_pix2coords_internal,    4},
  {"_rcosmo_pointInConvexPolygon",   (DL_FUNC) &_rcosmo_pointInConvexPolygon,   2},
  {"_rcosmo_pointInConvexPolygonHP", (DL_FUNC) &_rcosmo_pointInConvexPolygonHP, 4},
  {"_rcosmo_pointInDisc",            (DL_FUNC) &_rcosmo_pointInDisc,            2},
  {"_rcosmo_pointInDiscHP",          (DL_FUNC) &_rcosmo_pointInDiscHP,          4},
  {"_rcosmo_sph2car",                (DL_FUNC) &_rcosmo_sph2car,                1},
  {NULL, NULL, 0}
};

void R_init_rcosmo(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
