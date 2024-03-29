// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// covCMB_internal2
NumericVector covCMB_internal2(Rcpp::DataFrame cmbdf, NumericVector cos_breaks);
RcppExport SEXP _rcosmo_covCMB_internal2(SEXP cmbdfSEXP, SEXP cos_breaksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type cmbdf(cmbdfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cos_breaks(cos_breaksSEXP);
    rcpp_result_gen = Rcpp::wrap(covCMB_internal2(cmbdf, cos_breaks));
    return rcpp_result_gen;
END_RCPP
}
// covCMB_internal_var
NumericVector covCMB_internal_var(Rcpp::DataFrame cmbdf, NumericVector cos_breaks);
RcppExport SEXP _rcosmo_covCMB_internal_var(SEXP cmbdfSEXP, SEXP cos_breaksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type cmbdf(cmbdfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cos_breaks(cos_breaksSEXP);
    rcpp_result_gen = Rcpp::wrap(covCMB_internal_var(cmbdf, cos_breaks));
    return rcpp_result_gen;
END_RCPP
}
// mkpix2xyC
IntegerMatrix mkpix2xyC(int nside);
RcppExport SEXP _rcosmo_mkpix2xyC(SEXP nsideSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nside(nsideSEXP);
    rcpp_result_gen = Rcpp::wrap(mkpix2xyC(nside));
    return rcpp_result_gen;
END_RCPP
}
// nest2ring
IntegerVector nest2ring(int nside, IntegerVector pix);
RcppExport SEXP _rcosmo_nest2ring(SEXP nsideSEXP, SEXP pixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nside(nsideSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pix(pixSEXP);
    rcpp_result_gen = Rcpp::wrap(nest2ring(nside, pix));
    return rcpp_result_gen;
END_RCPP
}
// pix2coords_internal
NumericMatrix pix2coords_internal(int nside, bool nested, Rcpp::Nullable<Rcpp::IntegerVector> spix, bool cartesian);
RcppExport SEXP _rcosmo_pix2coords_internal(SEXP nsideSEXP, SEXP nestedSEXP, SEXP spixSEXP, SEXP cartesianSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nside(nsideSEXP);
    Rcpp::traits::input_parameter< bool >::type nested(nestedSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type spix(spixSEXP);
    Rcpp::traits::input_parameter< bool >::type cartesian(cartesianSEXP);
    rcpp_result_gen = Rcpp::wrap(pix2coords_internal(nside, nested, spix, cartesian));
    return rcpp_result_gen;
END_RCPP
}
// car2sph
DataFrame car2sph(DataFrame df);
RcppExport SEXP _rcosmo_car2sph(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(car2sph(df));
    return rcpp_result_gen;
END_RCPP
}
// sph2car
DataFrame sph2car(DataFrame df);
RcppExport SEXP _rcosmo_sph2car(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(sph2car(df));
    return rcpp_result_gen;
END_RCPP
}
// pointInConvexPolygonHP
LogicalVector pointInConvexPolygonHP(int nside, bool nested, DataFrame win, Rcpp::Nullable<Rcpp::IntegerVector> spix);
RcppExport SEXP _rcosmo_pointInConvexPolygonHP(SEXP nsideSEXP, SEXP nestedSEXP, SEXP winSEXP, SEXP spixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nside(nsideSEXP);
    Rcpp::traits::input_parameter< bool >::type nested(nestedSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type win(winSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type spix(spixSEXP);
    rcpp_result_gen = Rcpp::wrap(pointInConvexPolygonHP(nside, nested, win, spix));
    return rcpp_result_gen;
END_RCPP
}
// pointInDiscHP
LogicalVector pointInDiscHP(int nside, bool nested, DataFrame win, Rcpp::Nullable<Rcpp::IntegerVector> spix);
RcppExport SEXP _rcosmo_pointInDiscHP(SEXP nsideSEXP, SEXP nestedSEXP, SEXP winSEXP, SEXP spixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nside(nsideSEXP);
    Rcpp::traits::input_parameter< bool >::type nested(nestedSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type win(winSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type spix(spixSEXP);
    rcpp_result_gen = Rcpp::wrap(pointInDiscHP(nside, nested, win, spix));
    return rcpp_result_gen;
END_RCPP
}
// pointInConvexPolygon
LogicalVector pointInConvexPolygon(DataFrame df, DataFrame win);
RcppExport SEXP _rcosmo_pointInConvexPolygon(SEXP dfSEXP, SEXP winSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type win(winSEXP);
    rcpp_result_gen = Rcpp::wrap(pointInConvexPolygon(df, win));
    return rcpp_result_gen;
END_RCPP
}
// pointInDisc
LogicalVector pointInDisc(DataFrame df, DataFrame win);
RcppExport SEXP _rcosmo_pointInDisc(SEXP dfSEXP, SEXP winSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type win(winSEXP);
    rcpp_result_gen = Rcpp::wrap(pointInDisc(df, win));
    return rcpp_result_gen;
END_RCPP
}
// minDist_internal1
double minDist_internal1(Rcpp::DataFrame cmbdf, NumericVector point);
RcppExport SEXP _rcosmo_minDist_internal1(SEXP cmbdfSEXP, SEXP pointSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type cmbdf(cmbdfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type point(pointSEXP);
    rcpp_result_gen = Rcpp::wrap(minDist_internal1(cmbdf, point));
    return rcpp_result_gen;
END_RCPP
}
// minDist_internal2
double minDist_internal2(Rcpp::DataFrame cmbdf);
RcppExport SEXP _rcosmo_minDist_internal2(SEXP cmbdfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type cmbdf(cmbdfSEXP);
    rcpp_result_gen = Rcpp::wrap(minDist_internal2(cmbdf));
    return rcpp_result_gen;
END_RCPP
}
// maxDist_internal1
double maxDist_internal1(Rcpp::DataFrame cmbdf, NumericVector point);
RcppExport SEXP _rcosmo_maxDist_internal1(SEXP cmbdfSEXP, SEXP pointSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type cmbdf(cmbdfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type point(pointSEXP);
    rcpp_result_gen = Rcpp::wrap(maxDist_internal1(cmbdf, point));
    return rcpp_result_gen;
END_RCPP
}
// maxDist_internal2
double maxDist_internal2(Rcpp::DataFrame cmbdf);
RcppExport SEXP _rcosmo_maxDist_internal2(SEXP cmbdfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type cmbdf(cmbdfSEXP);
    rcpp_result_gen = Rcpp::wrap(maxDist_internal2(cmbdf));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rcosmo_covCMB_internal2", (DL_FUNC) &_rcosmo_covCMB_internal2, 2},
    {"_rcosmo_covCMB_internal_var", (DL_FUNC) &_rcosmo_covCMB_internal_var, 2},
    {"_rcosmo_mkpix2xyC", (DL_FUNC) &_rcosmo_mkpix2xyC, 1},
    {"_rcosmo_nest2ring", (DL_FUNC) &_rcosmo_nest2ring, 2},
    {"_rcosmo_pix2coords_internal", (DL_FUNC) &_rcosmo_pix2coords_internal, 4},
    {"_rcosmo_car2sph", (DL_FUNC) &_rcosmo_car2sph, 1},
    {"_rcosmo_sph2car", (DL_FUNC) &_rcosmo_sph2car, 1},
    {"_rcosmo_pointInConvexPolygonHP", (DL_FUNC) &_rcosmo_pointInConvexPolygonHP, 4},
    {"_rcosmo_pointInDiscHP", (DL_FUNC) &_rcosmo_pointInDiscHP, 4},
    {"_rcosmo_pointInConvexPolygon", (DL_FUNC) &_rcosmo_pointInConvexPolygon, 2},
    {"_rcosmo_pointInDisc", (DL_FUNC) &_rcosmo_pointInDisc, 2},
    {"_rcosmo_minDist_internal1", (DL_FUNC) &_rcosmo_minDist_internal1, 2},
    {"_rcosmo_minDist_internal2", (DL_FUNC) &_rcosmo_minDist_internal2, 1},
    {"_rcosmo_maxDist_internal1", (DL_FUNC) &_rcosmo_maxDist_internal1, 2},
    {"_rcosmo_maxDist_internal2", (DL_FUNC) &_rcosmo_maxDist_internal2, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_rcosmo(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
