#include "utils.hpp"
#pragma once

struct DemogParam {
  epp::map_of_matrix<double> basepop;
  double *    births;
  epp::map_of_cube<double> Sx;
  epp::map_of_cube<double> netmigr;
  epp::map_of_matrix<double> asfr;
  epp::map_of_matrix<double> srb;
  epp::map_of_matrix<double> birthslag;
  epp::map_of_matrix<double> cumsurv;
  epp::map_of_matrix<double> cumnetmigr;
  epp::map_of_cube<double> targetpop;
  epp::map_of_matrix<double> entrantpop;
  const bool  flag_popadjust;
  DemogParam(const SEXP& fp) :
    basepop   (REAL(get_value(fp, "basepop")), get_dim_2D(fp, "basepop")),
    births    (REAL(get_value(fp, "births"))),
    Sx        (REAL(get_value(fp, "Sx")), get_dim_3D(fp, "Sx")),
    netmigr   (REAL(get_value(fp, "netmigr")), get_dim_3D(fp, "netmigr")),
    asfr      (REAL(get_value(fp, "asfr")), get_dim_2D(fp, "asfr")),
    srb       (REAL(get_value(fp, "srb")), get_dim_2D(fp, "srb")),
    birthslag (REAL(get_value(fp, "birthslag")), get_dim_2D(fp, "birthslag")),
    cumsurv   (REAL(get_value(fp, "cumsurv")), get_dim_2D(fp, "cumsurv")),
    cumnetmigr(REAL(get_value(fp, "cumnetmigr")), get_dim_2D(fp, "cumnetmigr")),
    targetpop (REAL(get_value(fp, "targetpop")), get_dim_3D(fp, "targetpop")),
    entrantpop(REAL(get_value(fp, "entrantpop")), get_dim_2D(fp, "entrantpop")),
    flag_popadjust ((bool) *REAL(get_value(fp, "popadjust")))
    {}

};

struct NaturalHistoryParam {
  epp::map_of_matrix<double> artmx_timerr;
  epp::map_of_cube<double> cd4_initdist;
  epp::map_of_cube<double> cd4_prog;
  epp::map_of_cube<double> cd4_mort;
  epp::map_of_cube<double> frr_cd4;
  epp::map_of_4D<double> art_mort;
  epp::map_of_4D<double> frr_art;
  NaturalHistoryParam(const SEXP& fp) :
    artmx_timerr(REAL(get_value(fp, "artmx_timerr")), get_dim_2D(fp, "artmx_timerr")),
    cd4_initdist(REAL(get_value(fp, "cd4_initdist")), get_dim_3D(fp, "cd4_initdist")),
    cd4_prog    (REAL(get_value(fp, "cd4_prog")), get_dim_3D(fp, "cd4_prog")),
    cd4_mort    (REAL(get_value(fp, "cd4_mort")), get_dim_3D(fp, "cd4_mort")),
    frr_cd4     (REAL(get_value(fp, "frr_cd4")), get_dim_3D(fp, "frr_cd4")),
    art_mort    (REAL(get_value(fp, "art_mort")), get_dim_4D(fp, "art_mort")),
    frr_art     (REAL(get_value(fp, "frr_art")), get_dim_4D(fp, "frr_art"))
    {}
};

struct ArtData {
  epp::map_of_matrix<double>    art15plus_num;
  epp::map_of_matrix<double>    art15plus_isperc;
  const double * artcd4elig_idx;  // NOTE: 1-based indexing
  const double * specpop_percelig;
  const double * pw_artelig;
  const double   who34percelig;
  const double * art_dropout;
  const double * median_cd4init;
  const double * med_cd4init_cat;
  const double * med_cd4init_input;
  const int      art_alloc_method;
  const double   art_alloc_mxweight;
  const int      scale_cd4_mort;
  ArtData(const SEXP& fp) :
    art15plus_num      (REAL(get_value(fp, "art15plus_num")),
                             get_dim_2D(fp, "art15plus_num")),
    art15plus_isperc   (REAL(get_value(fp, "art15plus_isperc")),
                             get_dim_2D(fp, "art15plus_isperc")),
    artcd4elig_idx     (REAL(get_value(fp, "artcd4elig_idx"))),
    specpop_percelig   (REAL(get_value(fp, "specpop_percelig"))),
    pw_artelig         (REAL(get_value(fp, "pw_artelig"))),
    who34percelig      (*REAL(get_value(fp, "who34percelig"))),
    art_dropout        (REAL(get_value(fp, "art_dropout"))),
    median_cd4init     (REAL(get_value(fp, "median_cd4init"))),
    med_cd4init_cat    (REAL(get_value(fp, "med_cd4init_cat"))),
    med_cd4init_input  (REAL(get_value(fp, "med_cd4init_input"))),
    art_alloc_method   ((int) *REAL(get_value(fp, "art_alloc_method"))),
    art_alloc_mxweight (*REAL(get_value(fp, "art_alloc_mxweight"))),
    scale_cd4_mort     ((int) *REAL(get_value(fp, "scale_cd4_mort")))
    {}
};

struct RtrendParam {
  const double * proj_steps;
  const double * rw_start;
  const double * rw_trans;
  const double * rlogistic_steps;
  const double * rw_steps;
  const double * n_rw;
  const double * rw_dk;
  const double * rw_knots;
  const double * rw_idx;
  const double * n_param;
  const double * rw_transition;
  RtrendParam() {};
  void init_me(const SEXP& fp) {
    SEXP fp_rt      = get_value(fp, "rt");
    proj_steps      = REAL(get_value(fp_rt, "proj_steps"));
    rw_start        = REAL(get_value(fp_rt, "rw_start"));
    rw_trans        = REAL(get_value(fp_rt, "rw_trans"));
    rlogistic_steps = REAL(get_value(fp_rt, "rlogistic_steps"));
    rw_steps        = REAL(get_value(fp_rt, "rw_steps"));
    n_rw            = REAL(get_value(fp_rt, "n_rw"));
    rw_dk           = REAL(get_value(fp_rt, "rw_dk"));
    rw_knots        = REAL(get_value(fp_rt, "rw_knots"));
    rw_idx          = REAL(get_value(fp_rt, "rw_idx"));
    n_param         = REAL(get_value(fp_rt, "n_param"));
    rw_transition   = REAL(get_value(fp_rt, "rw_transition"));
  }
};

struct IncidenceParam {
  const int      eppmod; 
  const int      incidmod;
  epp::map_of_cube<double> incrr_age;
  epp::map_of_matrix<double> circ_prop;
  epp::map_of_cube<double> mixmat;
  epp::map_of_cube<double> db_rate;
  epp::map_of_cube<double> est_condom;
  epp::map_of_matrix<double> est_senesence;
  epp::map_of_matrix<double> est_pcr;
  epp::map_of_matrix<double> leading_ev;
  const double   balancing;
  const double   stage0_time;
  const double   relinfectART;
  const double * incrr_sex;
  const double * mf_transm_rr;
  const double * rel_vl;
  epp::map_of_matrix<double> fage;
  double         circ_incid_rr;
  const double * incidinput;
  int            incidpopage;
  double         tsEpidemicStart; //ts_epidemic_start;
  double       * proj_steps;
  double         iota;
  const double * logitiota;
  epp::map_of_vector<double> rvec;
  double         rw_start;
  RtrendParam    rt;
  IncidenceParam(const SEXP& fp) :
      eppmod        ((int) *REAL(get_value(fp, "eppmodInt"))),
      incidmod      ((int) *REAL(get_value(fp, "incidmodInt"))),
      incrr_age     ( REAL(get_value(fp, "incrr_age")), get_dim_3D(fp, "incrr_age")),
      circ_prop     ( REAL(get_value(fp, "circ_prop")), get_dim_2D(fp, "circ_prop")),
      mixmat        ( REAL(get_value(fp, "mixmat")), get_dim_3D(fp, "mixmat")),
      db_rate       ( REAL(get_value(fp, "db_rate")), get_dim_3D(fp, "db_rate")),
      est_condom    ( REAL(get_value(fp, "est_condom")), get_dim_3D(fp, "est_condom")),
      est_senesence ( REAL(get_value(fp, "est_senesence")), get_dim_2D(fp, "est_senesence")),
      est_pcr       ( REAL(get_value(fp, "est_pcr")), get_dim_2D(fp, "est_pcr")),
      leading_ev    ( REAL(get_value(fp, "leading_ev")), get_dim_2D(fp, "leading_ev")),
      balancing     (*REAL(get_value(fp, "balancing"))),
      stage0_time   (*REAL(get_value(fp, "stage0_time"))),
      relinfectART  (*REAL(get_value(fp, "relinfectART"))),
      incrr_sex     ( REAL(get_value(fp, "incrr_sex"))),
      mf_transm_rr  ( REAL(get_value(fp, "mf_transm_rr"))),
      rel_vl        ( REAL(get_value(fp, "rel_vl"))),
      fage          ( REAL(get_value(fp, "fage")), get_dim_2D(fp, "fage")),
      circ_incid_rr (*REAL(get_value(fp, "circ_incid_rr"))),
      rvec          ( REAL(get_value(fp, "rvec"))         , get_dim_1D(fp, "rvec"))
    {
      if (eppmod == 2) {  // direct incidence input
        incidinput     =        REAL(get_value(fp, "incidinput"));
        incidpopage    = (int) *REAL(get_value(fp, "incidpopage"));
      }
      if (eppmod != 2) {  // != direct incidence input
        tsEpidemicStart  = *REAL(get_value(fp, "tsEpidemicStart"));
        iota             = *REAL(get_value(fp, "iota"));
        logitiota        =  REAL(get_value(fp, "logitiota"));
        proj_steps       =  REAL(get_value(fp, "proj_steps"));
      }
      if (eppmod == 0)  // rhybrid
        rw_start         = *REAL(get_value(fp, "rw_start"));
      if (eppmod == 1)  // rtrend
        rt.init_me(fp);
    }
};

struct PaediatricHivParam {
  const double * verttrans_lag;
  const double * paedsurv_lag;
  const double   netmighivsurv;
  const double   netmig_hivprob;
  epp::map_of_matrix<double>    entrantprev;
  epp::map_of_matrix<double>    entrantartcov;
  epp::map_of_cube<double>    paedsurv_cd4dist;
  epp::map_of_4D<double>    paedsurv_artcd4dist;
  PaediatricHivParam(const SEXP& fp) :
      verttrans_lag  (REAL(get_value(fp, "verttrans_lag"))),
      paedsurv_lag   (REAL(get_value(fp, "paedsurv_lag"))),
      netmighivsurv  (*REAL(get_value(fp, "netmighivsurv"))),
      netmig_hivprob (*REAL(get_value(fp, "netmig_hivprob"))),
      entrantprev         (REAL(get_value(fp, "entrantprev")),
                           get_dim_2D(fp, "entrantprev")),
      entrantartcov       (REAL(get_value(fp, "entrantartcov")),
                           get_dim_2D(fp, "entrantartcov")),
      paedsurv_cd4dist    (REAL(get_value(fp, "paedsurv_cd4dist")),
                           get_dim_3D(fp, "paedsurv_cd4dist")),
      paedsurv_artcd4dist (REAL(get_value(fp, "paedsurv_artcd4dist")),
                           get_dim_4D(fp, "paedsurv_artcd4dist"))
    {}
};

struct AncParam {
  const double * ancsitedata;
  const double * ancrt;
  const double * ancbias;
  const double * v_infl;
  const double * ancrtcens_vinfl;
  const double * ancrtsite_beta;
  const double * log_frr_adjust;
  AncParam(const SEXP& fp) {
      if (has_value(fp, "ancsitedata")) {
        ancsitedata      = REAL(get_value(fp, "ancsitedata"));
        ancrt            = REAL(get_value(fp, "ancrtInt"));
      }
      if (has_value(fp, "ancbias")) {// double
        ancbias          = REAL(get_value(fp, "ancbias"));      
        v_infl           = REAL(get_value(fp, "v_infl"));      
        ancrtcens_vinfl  = REAL(get_value(fp, "ancrtcens_vinfl"));      
        ancrtsite_beta   = REAL(get_value(fp, "ancrtsite_beta"));
        log_frr_adjust   = REAL(get_value(fp, "log_frr_adjust"));
      }    
    }
};

struct Parameters {
  DemogParam          dm;
  NaturalHistoryParam nh;
  ArtData             ad;
  IncidenceParam      ic;
  PaediatricHivParam  ph;
  Parameters(const SEXP& fp) :
    dm(fp), // DemogParam 
    nh(fp), // NaturalHistoryParam
    ad(fp), // ArtData
    ic(fp), // IncidenceParam
    ph(fp) // PaediatricHivParam
    {}
};

// Master parameters class
struct StateSpace {
  const int    SIM_YEARS;
  SEXP         fp_ss;
  int          year = 0; // simulation year TODO: should move out of StateSpace
  const int    MODEL;
  const bool   MIX;
  const int    proj_start;
  const int    PROJ_YEARS;
  const int    AGE_START;
  const int    steps_per_year;
  const int    NG;
  const int    pDS;
  const int    M;
  const int    F;
  const int    P; // positive
  const int    N; // negative
  const int    pAG;
  const double ag_rate;
  const int    hAG;
  const int    hDS;
  const int    hTS;
  const double DT;
  const int    pDB;
  const int    hDB;
  const int    n_steps;
  const int    tARTstart;
  epp::map_of_vector<double> p_fert_;
  epp::map_of_vector<double> p_age15to49_;
  epp::map_of_vector<double> p_age15plus_;
  epp::map_of_vector<double> h_ag_span;
  epp::map_of_vector<double> ag_;
  epp::map_of_vector<double> agfirst_;
  epp::map_of_vector<double> aglast_;
  epp::map_of_vector<double> h_fert_;
  epp::map_of_vector<double> h_age15to49_;
  epp::map_of_vector<double> h_age15plus_;
  const int pAG_FERT, hAG_FERT, pAG_1549, hAG_1549, pAG_15plus, hAG_15plus;
  // 
  int art_dim() {return hTS * hDS * hAG * NG;};
  int hiv_dim() {return hDS * hAG * NG;};
  int pop_dim() {return pDS * pAG * NG;};
  StateSpace(const SEXP& fp) :
    SIM_YEARS      ((int) *REAL(get_value(fp, "SIM_YEARS"))),
    fp_ss          (            get_value(fp, "ss")),
    MODEL          ((int) *REAL(get_value(fp_ss, "MODEL"))),
    MIX            ((bool)*REAL(get_value(fp_ss, "MIX"))),
    proj_start     ((int) *REAL(get_value(fp_ss, "proj_start"))),
    PROJ_YEARS     ((int) *REAL(get_value(fp_ss, "PROJ_YEARS"))),
    AGE_START      ((int) *REAL(get_value(fp_ss, "AGE_START"))),
    steps_per_year ((int) *REAL(get_value(fp_ss, "hiv_steps_per_year"))),
    NG             ((int) *REAL(get_value(fp_ss, "NG"))),
    pDS            ((int) *REAL(get_value(fp_ss, "pDS"))),
    M              ((int) *REAL(get_value(fp_ss, "m_idx")) - 1),
    F              ((int) *REAL(get_value(fp_ss, "f_idx")) - 1),
    P              ((int) *REAL(get_value(fp_ss, "hivp_idx")) - 1),
    N              ((int) *REAL(get_value(fp_ss, "hivn_idx")) - 1),
    pAG            ((int) *REAL(get_value(fp_ss, "pAG"))),
    ag_rate        (      *REAL(get_value(fp_ss, "ag_rate"))),
    hAG            ((int) *REAL(get_value(fp_ss, "hAG"))),
    hDS            ((int) *REAL(get_value(fp_ss, "hDS"))),
    hTS            ((int) *REAL(get_value(fp_ss, "hTS"))),
    DT             (      *REAL(get_value(fp_ss, "DT"))),
    pDB            ((int) *REAL(get_value(fp_ss, "pDB"))),
    hDB            (pDB),
    n_steps        ((int) *REAL(get_value(fp, "n_steps"))),
    tARTstart      ((int) *REAL(get_value(fp, "tARTstart"))),
    p_fert_        (REAL(get_value(fp_ss, "p_fert_idx")), get_dim_1D(fp_ss, "p_fert_idx")),
    p_age15to49_   (REAL(get_value(fp_ss, "p_age15to49_idx")), get_dim_1D(fp_ss, "p_age15to49_idx")),
    p_age15plus_   (REAL(get_value(fp_ss, "p_age15plus_idx")), get_dim_1D(fp_ss, "p_age15plus_idx")),
    h_ag_span      (REAL(get_value(fp_ss, "h_ag_span")), get_dim_1D(fp_ss, "h_ag_span")),
    ag_            (REAL(get_value(fp_ss, "ag_idx")), get_dim_1D(fp_ss, "ag_idx")),
    agfirst_       (REAL(get_value(fp_ss, "agfirst_idx")), get_dim_1D(fp_ss, "agfirst_idx")),
    aglast_        (REAL(get_value(fp_ss, "aglast_idx")), get_dim_1D(fp_ss, "aglast_idx")),
    h_fert_        (REAL(get_value(fp_ss, "h_fert_idx")), get_dim_1D(fp_ss, "h_fert_idx")),
    h_age15to49_   (REAL(get_value(fp_ss, "h_age15to49_idx")), get_dim_1D(fp_ss, "h_age15to49_idx")),
    h_age15plus_   (REAL(get_value(fp_ss, "h_age15plus_idx")), get_dim_1D(fp_ss, "h_age15plus_idx")),
    pAG_FERT       ((p_fert_(0)      - 1) + p_fert_.size()),
    hAG_FERT       ((h_fert_(0)      - 1) + h_fert_.size()),
    pAG_1549       ((p_age15to49_(0) - 1) + p_age15to49_.size()),
    hAG_1549       ((h_age15to49_(0) - 1) + h_age15to49_.size()),
    pAG_15plus     ((p_age15plus_(0) - 1) + p_age15plus_.size()),
    hAG_15plus     ((h_age15plus_(0) - 1) + h_age15plus_.size())
    {}
};

template<typename T>
struct Views { // two years view of outputs
  const int         N_pop;
  const int         N_hiv;
  const int         N_art;
  const std::array<long, 3> pop_shape;
  const std::array<long, 3> hiv_shape;
  const std::array<long, 4> art_shape;
  epp::map_of_cube<T> now_pop;
  epp::map_of_cube<T> now_hiv;
  epp::map_of_4D<T> now_art;
  epp::map_of_cube<T> pre_pop; // last year view
  epp::map_of_cube<T> pre_hiv; // last year view
  epp::map_of_4D<T> pre_art; // last year view
  Views(T * pop_start,
        T * hiv_start,
        T * art_start,
        StateSpace* const s, int year = 1) :
    N_pop     (s->NG * s->pAG * s->pDS),
    N_hiv     (s->NG * s->hAG * s->hDS),
    N_art     (s->NG * s->hAG * s->hDS * s->hTS),
    pop_shape ({ s->pAG, s->NG, s->pDS }),
    hiv_shape ({         s->hDS, s->hAG, s->NG}),
    art_shape ({ s->hTS, s->hDS, s->hAG, s->NG}),
    now_pop   (pop_start + year * N_pop, pop_shape),
    now_hiv   (hiv_start + year * N_hiv, hiv_shape),
    now_art   (art_start + year * N_art, art_shape),
    pre_pop   (pop_start + (year - 1) * N_pop, pop_shape),
    pre_hiv   (hiv_start + (year - 1) * N_hiv, hiv_shape),
    pre_art   (art_start + (year - 1) * N_art, art_shape)
  {}
};