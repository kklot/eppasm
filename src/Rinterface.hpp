#pragma once
#include "fpClass.hpp"
#include "utils.hpp"
#include "model.hpp"

template<typename T>
struct outputSEXP {
  SEXP 
    pop,
    hivpop,
    artpop,
    data_db,
    hivpopdb,
    artpopdb,
    prev15to49,
    incid15to49,
    entrantprev,
    pregprevlag,
    pregprev,
    inci15to49_ts,
    prev15to49_ts,
    rvec,
    infections,
    hivdeaths,
    natdeaths,
    popadjust,
    stage0;
  int np = 0;
  outputSEXP(const StateSpace& s, const Model<T> * m) {
    // Allocate memory
    int size_art = s.hTS * s.hDS * s.hAG * s.NG * s.PROJ_YEARS,
        size_hiv =         s.hDS * s.hAG * s.NG * s.PROJ_YEARS, 
        size_pop = s.pAG * s.NG * s.pDS * s.PROJ_YEARS, 
        size_pop_vg = s.pDB * s.NG * s.pDS * s.PROJ_YEARS, 
        size_age_sex_year = s.pAG * s.NG * s.PROJ_YEARS;

    artpop        = PROTECT(NEW_NUMERIC(size_art)); ++np;
    artpopdb      = PROTECT(NEW_NUMERIC(size_art)); ++np;
    hivpop        = PROTECT(NEW_NUMERIC(size_hiv)); ++np;
    hivpopdb      = PROTECT(NEW_NUMERIC(size_hiv)); ++np;
    pop           = PROTECT(NEW_NUMERIC(size_pop)); ++np;
    prev15to49    = PROTECT(NEW_NUMERIC(s.PROJ_YEARS)); ++np;
    incid15to49   = PROTECT(NEW_NUMERIC(s.PROJ_YEARS)); ++np;
    entrantprev   = PROTECT(NEW_NUMERIC(s.PROJ_YEARS)); ++np;
    pregprevlag   = PROTECT(NEW_NUMERIC(s.PROJ_YEARS)); ++np;
    pregprev      = PROTECT(NEW_NUMERIC(s.PROJ_YEARS)); ++np;
    inci15to49_ts = PROTECT(NEW_NUMERIC(s.n_steps)); ++np;
    prev15to49_ts = PROTECT(NEW_NUMERIC(s.n_steps)); ++np;
    rvec          = PROTECT(NEW_NUMERIC(s.n_steps)); ++np;
    infections    = PROTECT(NEW_NUMERIC(size_age_sex_year)); ++np;
    hivdeaths     = PROTECT(NEW_NUMERIC(size_age_sex_year)); ++np;
    natdeaths     = PROTECT(NEW_NUMERIC(size_age_sex_year)); ++np;
    popadjust     = PROTECT(NEW_NUMERIC(size_age_sex_year)); ++np;    
    data_db       = PROTECT(NEW_NUMERIC(size_pop_vg)); ++np;
    stage0        = PROTECT(NEW_NUMERIC(s.hAG * s.NG * s.PROJ_YEARS)); ++np;

    // Copying this because we want the model to be independent from R, otherwise
    // we can map and work directly on the memory above
    memcpy(REAL(stage0)        , m->hivpop.stage0.data()        , sizeof(T) * s.hAG * s.NG * s.PROJ_YEARS);
    memcpy(REAL(hivpopdb)      , m->hivpop.data_db.data()       , sizeof(T) * size_hiv);
    memcpy(REAL(hivpop)        , m->hivpop.data.data()          , sizeof(T) * size_hiv);
    memcpy(REAL(artpopdb)      , m->artpop.data_db.data()       , sizeof(T) * size_art);
    memcpy(REAL(artpop)        , m->artpop.data.data()          , sizeof(T) * size_art);
    memcpy(REAL(prev15to49)    , m->pop.prev15to49.data()       , sizeof(T) * s.PROJ_YEARS);
    memcpy(REAL(incid15to49)   , m->pop.incid15to49.data()      , sizeof(T) * s.PROJ_YEARS);
    memcpy(REAL(entrantprev)   , m->pop.entrantprev.data()      , sizeof(T) * s.PROJ_YEARS);
    memcpy(REAL(pregprevlag)   , m->pop.pregprevlag.data()      , sizeof(T) * s.PROJ_YEARS);
    memcpy(REAL(pregprev)      , m->pop.pregprev.data()         , sizeof(T) * s.PROJ_YEARS);
    memcpy(REAL(inci15to49_ts) , m->pop.incrate15to49_ts.data() , sizeof(T) * s.n_steps);
    memcpy(REAL(prev15to49_ts) , m->pop.prev15to49_ts.data()    , sizeof(T) * s.n_steps);
    memcpy(REAL(rvec)          , m->pop.rvec.data()             , sizeof(T) * s.n_steps);
    memcpy(REAL(infections)    , m->pop.infections.data()       , sizeof(T) * size_age_sex_year);
    memcpy(REAL(hivdeaths)     , m->pop.hivdeaths.data()        , sizeof(T) * size_age_sex_year);
    memcpy(REAL(natdeaths)     , m->pop.natdeaths.data()        , sizeof(T) * size_age_sex_year);
    memcpy(REAL(popadjust)     , m->pop.popadjust.data()        , sizeof(T) * size_age_sex_year);
    memcpy(REAL(pop)           , m->pop.data.data()             , sizeof(T) * size_pop);
    memcpy(REAL(data_db)       , m->pop.data_db.data()          , sizeof(T) * size_pop_vg);
    
    // set dimensions for R
    SEXP pop_sexp_dim = PROTECT(NEW_INTEGER(4)); ++np;
    INTEGER(pop_sexp_dim)[0] = s.pAG;
    INTEGER(pop_sexp_dim)[1] = s.NG;
    INTEGER(pop_sexp_dim)[2] = s.pDS;
    INTEGER(pop_sexp_dim)[3] = s.PROJ_YEARS;
    SET_DIM(pop, pop_sexp_dim);

    if (s.MODEL!=0) {
      SEXP age_sex_year_dim = PROTECT(NEW_INTEGER(3)); ++np;
      INTEGER(age_sex_year_dim)[0] = s.pAG;
      INTEGER(age_sex_year_dim)[1] = s.NG;
      INTEGER(age_sex_year_dim)[2] = s.PROJ_YEARS;
      SET_DIM(infections, age_sex_year_dim);
      SET_DIM(hivdeaths, age_sex_year_dim);
      SET_DIM(natdeaths, age_sex_year_dim);
      SET_DIM(popadjust, age_sex_year_dim);

      SEXP hiv_sexp_dim = PROTECT(NEW_INTEGER(4)); ++np;
      INTEGER(hiv_sexp_dim)[0] = s.hDS;
      INTEGER(hiv_sexp_dim)[1] = s.hAG;
      INTEGER(hiv_sexp_dim)[2] = s.NG;
      INTEGER(hiv_sexp_dim)[3] = s.PROJ_YEARS;
      SET_DIM(hivpop, hiv_sexp_dim);

      SEXP art_sexp_dim = PROTECT(NEW_INTEGER(5)); ++np;
      INTEGER(art_sexp_dim)[0] = s.hTS;
      INTEGER(art_sexp_dim)[1] = s.hDS;
      INTEGER(art_sexp_dim)[2] = s.hAG;
      INTEGER(art_sexp_dim)[3] = s.NG;
      INTEGER(art_sexp_dim)[4] = s.PROJ_YEARS;
      SET_DIM(artpop, art_sexp_dim);

      SET_ATTR(pop, Rf_install("infections"), infections);
      SET_ATTR(pop, Rf_install("hivdeaths"), hivdeaths);
      SET_ATTR(pop, Rf_install("natdeaths"), natdeaths);
      SET_ATTR(pop, Rf_install("popadjust"), popadjust);
      SET_ATTR(pop, Rf_install("pregprevlag"), pregprevlag);
      SET_ATTR(pop, Rf_install("pregprev"), pregprev);
      SET_ATTR(pop, Rf_install("rvec_ts"), rvec);
      SET_ATTR(pop, Rf_install("prev15to49"), prev15to49);
      SET_ATTR(pop, Rf_install("incid15to49"), incid15to49);
      SET_ATTR(pop, Rf_install("entrantprev"), entrantprev);
      SET_ATTR(pop, Rf_install("incrate15to49_ts"), inci15to49_ts);
      SET_ATTR(pop, Rf_install("prev15to49_ts_sexp"), prev15to49_ts);
      SET_ATTR(pop, Rf_install("artpop"), artpop);
      SET_ATTR(pop, Rf_install("hivpop"), hivpop);
      SET_CLASS(pop, Rf_mkString("spec"));

      if (s.MODEL == 2) {
        SEXP pop_db_sexp_dim = PROTECT(NEW_INTEGER(4)); ++np;
        INTEGER(pop_db_sexp_dim)[0] = s.pDB;
        INTEGER(pop_db_sexp_dim)[1] = s.NG;
        INTEGER(pop_db_sexp_dim)[2] = s.pDS;
        INTEGER(pop_db_sexp_dim)[3] = s.PROJ_YEARS;
        SET_DIM(data_db, pop_db_sexp_dim);

        SET_ATTR(pop, Rf_install("vpop"), data_db);
        
        SET_DIM(hivpopdb, hiv_sexp_dim);
        SET_DIM(artpopdb, art_sexp_dim);
        SET_ATTR(pop, Rf_install("vpopart"), artpopdb);
        SET_ATTR(pop, Rf_install("vpophiv"), hivpopdb);

        SEXP stage0_dim = PROTECT(NEW_INTEGER(3)); ++np;
        INTEGER(stage0_dim)[0] = s.hAG;
        INTEGER(stage0_dim)[1] = s.NG;
        INTEGER(stage0_dim)[2] = s.PROJ_YEARS;
        SET_DIM(stage0, stage0_dim);
        SET_ATTR(pop, Rf_install("stage0"), stage0);
      }
    }
    else 
      SET_CLASS(pop, Rf_mkString("dempp"));
    UNPROTECT(np); // pop 
  }
};