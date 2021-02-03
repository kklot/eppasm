#include "fpClass.hpp"

DemogParam::DemogParam(const SEXP& fp) :
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

NaturalHistoryParam::NaturalHistoryParam(const SEXP& fp) :
  artmx_timerr(REAL(get_value(fp, "artmx_timerr")), get_dim_2D(fp, "artmx_timerr")),
  cd4_initdist(REAL(get_value(fp, "cd4_initdist")), get_dim_3D(fp, "cd4_initdist")),
  cd4_prog    (REAL(get_value(fp, "cd4_prog")), get_dim_3D(fp, "cd4_prog")),
  cd4_mort    (REAL(get_value(fp, "cd4_mort")), get_dim_3D(fp, "cd4_mort")),
  frr_cd4     (REAL(get_value(fp, "frr_cd4")), get_dim_3D(fp, "frr_cd4")),
  art_mort    (REAL(get_value(fp, "art_mort")), get_dim_4D(fp, "art_mort")),
  frr_art     (REAL(get_value(fp, "frr_art")), get_dim_4D(fp, "frr_art"))
  {}

ArtData::ArtData(const SEXP& fp) :
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

void RtrendParam::init_me(const SEXP& fp) {
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

IncidenceParam::IncidenceParam(const SEXP& fp) :
    eppmod        ((int) *REAL(get_value(fp, "eppmodInt"))),
    incidmod      ((int) *REAL(get_value(fp, "incidmodInt"))),
    incrr_age     ( REAL(get_value(fp, "incrr_age")), get_dim_3D(fp, "incrr_age")),
    circ_prop     ( REAL(get_value(fp, "circ_prop")), get_dim_2D(fp, "circ_prop")),
    mixmat        ( REAL(get_value(fp, "mixmat")), get_dim_3D(fp, "mixmat")),
    db_rate       ( REAL(get_value(fp, "db_rate")), get_dim_3D(fp, "db_rate")),
    est_condom    ( REAL(get_value(fp, "est_condom")), get_dim_3D(fp, "est_condom")),
    est_senesence ( REAL(get_value(fp, "est_senesence")), get_dim_2D(fp, "est_senesence")),
    est_pcr       ( REAL(get_value(fp, "est_pcr")), get_dim_2D(fp, "est_pcr")),
    balancing     (*REAL(get_value(fp, "balancing"))),
    relinfectART  (*REAL(get_value(fp, "relinfectART"))),
    incrr_sex     ( REAL(get_value(fp, "incrr_sex"))),
    mf_transm_rr  ( REAL(get_value(fp, "mf_transm_rr"))),
    rel_vl  		  ( REAL(get_value(fp, "rel_vl"))),
    fage          ( REAL(get_value(fp, "fage")), get_dim_2D(fp, "fage")),
    circ_incid_rr (*REAL(get_value(fp, "circ_incid_rr")))
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
    if (has_value(fp, "rvec"))
      rvec             = REAL(get_value(fp, "rvec"));
    if (eppmod == 0)  // rhybrid
      rw_start         = *REAL(get_value(fp, "rw_start"));
    if (eppmod == 1)  // rtrend
      rt.init_me(fp);
  }

PaediatricHivParam::PaediatricHivParam(const SEXP& fp) :
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

AncParam::AncParam(const SEXP& fp) {
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

Parameters::Parameters(const SEXP& fp) :
  dm(fp), // DemogParam 
  nh(fp), // NaturalHistoryParam
  ad(fp), // ArtData
  ic(fp), // IncidenceParam
  ph(fp) // PaediatricHivParam
  {}

StateSpace::StateSpace(const SEXP& fp) :
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
  pAG_FERT       ((p_fert_[0]      - 1) + p_fert_.num_elements()),
  hAG_FERT       ((h_fert_[0]      - 1) + h_fert_.num_elements()),
  pAG_1549       ((p_age15to49_[0] - 1) + p_age15to49_.num_elements()),
  hAG_1549       ((h_age15to49_[0] - 1) + h_age15to49_.num_elements()),
  pAG_15plus     ((p_age15plus_[0] - 1) + p_age15plus_.num_elements()),
  hAG_15plus     ((h_age15plus_[0] - 1) + h_age15plus_.num_elements())
  {}
