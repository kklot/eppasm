#' @export
simmod.specfp <- function(fp) {
	fp  <- fp_fill_missing(fp)

  MODEL   = ifelse(is.null(fp$ss$MODEL), 1, fp$ss$MODEL)
  MIX     = ifelse(is.null(fp$ss$MIX), FALSE, fp$ss$MIX)
  VERSION = ifelse(is.null(fp$VERSION), 'R', fp$VERSION)
  
  if (VERSION != "R") {
    if (VERSION=="K") { # C++ classes
      fp  <- prepare_fp_for_Cpp(fp)
      mod <- .Call(eppasmOOpp, fp)
      return(.asList(mod))
    } 
    else { # keep this for tests
      mod <- .Call(eppasmC, fp)
      class(mod) <- "spec"
      return(.asList(mod))
    }
  }
  pop     <- popEPP$new(fp, MODEL, VERSION, MIX)
  hivpop  <- hivEPP$new(fp, MODEL)
  artpop  <- artEPP$new(fp, MODEL)
  for (i in 2:fp$SIM_YEARS) {
    pop$update_year <- hivpop$year <- artpop$year <- i
    epp_aging(pop, hivpop, artpop)
    epp_death(pop, hivpop, artpop)
    epp_migration(pop, hivpop, artpop)
    pop$update_fertile()
    if (MODEL != 0) { # Disease model simulation: events at dt timestep
      epp_disease_model(pop, hivpop, artpop)
      if (fp$eppmod == "directincid") ## Direct incidence input model
        pop$epp_disease_model_direct(hivpop, artpop)
    }
    if (fp$popadjust) # match target pop
      epp_adjustpop(pop, hivpop, artpop)
    if (MODEL != 0) {
      pop$cal_prev_pregant(hivpop, artpop) # prevalence among pregnant women
      pop$save_prev_n_inc() # save prevalence and incidence 15 to 49
    }
  }
  if (MODEL != 0) {
    pop$hivpop <- hivpop$data
    pop$artpop <- artpop$data
    if (MODEL==2) {
      pop$vpop    <- pop$VIRGIN$data
      pop$vpopart <- artpop$data_db
      pop$vpophiv <- hivpop$data_db
      pop$stage0  <- hivpop$stage0
    }
    class(pop) <- "spec"
  }
  else 
    class(pop) <- "dempp"
  return(pop)
}
