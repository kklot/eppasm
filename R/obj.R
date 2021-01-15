#' R6 class for a eppasm object - dev mode
#' 
#' Holding data, fitting methods,... 
#' 
#' @details
#' Read data, fit model variations, and visualize...
#' 
eppasm <- R6::R6Class('eppasm', lock_objects=FALSE, portable=FALSE, 
    public = list(
        data =  list(
            #' @field prev_15to49_nat country prevalence data
            prev_15to49_nat = NULL, 
            #' @field prev_agesex_nat country prevalence data
            prev_agesex_nat = NULL, 
            #' @field ancsitedata country prevalence data
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
            est_senesence = NULL
        ),
        #' @details
        #' Read default data, can be overrided by providing args to fits(...)
        #' 
        #' @param old_fits Reinitialize but copy old fitted objects, this is due
        #' to R6 once created reloading package does not update methods
        #' @return nothing
         
        initialize = function(old_fits) {
            if (!missing(old_fits)) 
                fits <<- old_fits
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
                obj           = data$inputs_nat[[country]],
                eppmod        = "rhybrid", 
                algorithm     = 'optim',
                fitincrr      = TRUE,
                rw_start      = 2005,
                version       = 'C'
            )
            attr(opts$obj, 'specfp')$relinfectART <- 1-0.7
            dbpar <- list(
                # sexual parameters
                db_rate       = data$est_db_rate[[country]],
                mixmat        = data$est_mixmat[[country]],
                est_senesence = data$est_senesence[[country]],
                est_pcr       = data$est_pcr[[country]],
                est_condom    = data$est_condom[[country]]
            ) # ignore in old code
            if (any(unlist(lapply(dbpar, is.null))))
                askYesNo('some sexual data is missing and filled by dummy value, continue?')
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
        },
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
            data$inputs_nat    <<- read_ext('inputs_nat.rds')
            # sexual parameters
            data$est_db_rate   <<- read_ext('est_db_fixed_first_year.rds')
            data$est_mixmat    <<- read_ext('est_mixmat_log_log_scaled.rds')
            data$est_pcr       <<- read_ext('est_pcr.rds')
            data$est_senesence <<- read_ext('est_senesence.rds')
            data$est_condom    <<- read_ext('est_condom.rds')
        },
        read_data = function() {
            data$prev_15to49_nat <<- read_ext('prev_15to49_nat.csv')
            data$prev_agesex_nat <<- read_ext('prev_agesex_nat.csv')
            data$ancsitedata     <<- read_ext('ancsitedata.csv')
        },
    ),
    active = list()
)