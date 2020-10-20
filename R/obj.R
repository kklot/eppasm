#' R6 class for a eppasm object - dev mode
#' 
#' Holding data, fitting methods,... 
#' 
#' @details
#' Read data, fit model variations, and visualize...
#' 
eppasm <- R6::R6Class('eppasm', lock_objects=FALSE, portable=FALSE, 
    public = list(
        #' @field data country prevalence data
        prev_15to49_nat = NULL, 
        #' @field data country prevalence data
        prev_agesex_nat = NULL, 
        #' @field data country prevalence data
        ancsitedata = NULL, 
        #' @field inputs_nat country prevalence data.
        inputs_nat = NULL, 
        #' @field est_db_rate input sexual debut rate.
        est_db_rate = NULL, 
        #' @field est_mixmat input sexual mixing matrices.
        est_mixmat = NULL, 
        #' @field est_pcr input partner change rate.
        est_pcr = NULL, 
        #' @field est_senesence input sexual senesence.
        est_senesence = NULL,

        #' @details
        #' Read default data, can be overrided by providing args to fits(...)
        #' 
        #' @param inputs no need for now
        #' @return nothing
         
        initialize = function(inputs) {
            private$read_input()
            private$read_data()
        }, 

        #' @field fits list of fitted object by name give in method \code{fit}
        fits = NULL,

        #' @details
        #' Call fitting
        #' 
        #' @param save_as name to save in \code{fits} member (list)
        #' @param country ISO code to choose country
        #' @param ... extra args to \code{\link{fitmod}}
        #' @return fit object save to `fits`
        #' @examples
        #' obj = eppasm$new()
        #' # run the fit
        #' obj$fit(save_as = 'default', country = 'MW')
        #' # access the fit
        #' obj$fits('default')

        fit = function(save_as, country, ...) {
            user <- list(...)
            opts <- list(
                # current model / remove eventually
                obj           = inputs_nat[[country]],
                eppmod        = "rhybrid", 
                algorithm     = 'optim',
                fitincrr      = TRUE,
                rw_start      = 2005,
                version       = 'C'
            )
            attr(opts$obj, 'specfp')$relinfectART <- 1-0.7
            dbpar <- list(
                # sexual parameters
                db_rate       = est_db_rate[[country]],
                mixmat        = est_mixmat[[country]],
                est_senesence = est_senesence[[country]], 
                est_pcr       = est_pcr[[country]],
                est_condom    = est_condom
            ) # ignore in old code
            # update options
            dbpar <- modifyList(dbpar, user)
            opts <- modifyList(opts, dbpar)
            fits[[save_as]] <<- do.call('fitmod', opts)
        },

        #' @details
        #' Save fit objects to files
        #' 
        #' @param name name of the fit in fits
        #' @param with_data save the default data?
        #' @return nothing
        backup = function(name, with_data=FALSE) {
            if (with_data)
                saveRDS(self, name)
            else
                saveRDS(fits, name)
        },

        #' @details
        #' simulate the model from the fits, use to e.g. debug as it changes the
        #' version to R not the C++ used to fit
        #' 
        #' @param name name of the fits
        #' @param version default to R, not C++ ("C" or "K")
        #' @return fitmod

        simulate = function(name, version="R") {
            with(fits[[name]], {
                fp <- modifyList(fp, list(VERSION=version))
                simmod(update(fp, list=fnCreateParam(par, fp)))
            })
        print = function(...) {
            cat("EPPASM object\ntype help(eppasm) for details...\n")
        }
    ),
    private = list(
        read_ext = function(name,...) {
            pth = system.file("extdata", name, package="eppasm")
            switch(tools::file_ext(pth),
                rds = readRDS(pth,...),
                csv = read.csv(pth,...)
            )
        },
        read_input = function(x=0) {
            # eppasm_inputs
            inputs_nat    <<- read_ext('inputs_nat.rds')
            # sexual parameters
            est_db_rate   <<- read_ext('est_db_fixed_first_year.rds')
            est_mixmat    <<- read_ext('est_mixmat_log_log_scaled.rds')
            est_pcr       <<- read_ext('est_pcr.rds')
            est_senesence <<- read_ext('est_senesence.rds')
            est_condom    <<- read_ext('est_condom_logistic_malawi.rds')
        },
        read_data = function() {
            prev_15to49_nat <<- read_ext('prev_15to49_nat.csv')
            prev_agesex_nat <<- read_ext('prev_agesex_nat.csv')
            ancsitedata     <<- read_ext('ancsitedata.csv')
        }
    ),
    active = list()
)