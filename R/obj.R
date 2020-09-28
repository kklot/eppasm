#' R6 class for a eppasm object 
#' 
#' Holding data, fitting methods,...
#' 
#' @details
#' Read data, fit model variations, and visualize...
#' 
eppasm <- R6::R6Class('eppasm', lock_objects=FALSE, portable=FALSE, 
    public = list(
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
        #' Read default data
        #' 
        #' @param inputs no need for now
        #' @return nothing
         
        initialize = function(inputs) {
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

        fit = function(save_as, country, mixing=FALSE, ...) {
            # current model
            opts <- list(
              obj           = inputs_nat[[country]],
              eppmod        = "rhybrid", 
              algorithm     = 'optim',
              fitincrr      = TRUE,
              rw_start      = 2005,
              version       = 'C'
            )
            attr(opts$obj, 'specfp')$relinfectART <- 1-0.7
            fits[[save_as]] <<- do.call('fitmod', opts)
        }
    ),
    private = list(
        read_ext = function(name) {
            pth = system.file("extdata", name, package="eppasm")
            readRDS(pth)
        },        
        read_data = function(x=0) {
            # eppasm_inputs
            inputs_nat <<- read_ext('inputs_nat.rds')
            # sexual parameters
            est_db_rate <<- read_ext('est_db_fixed_first_year.rds')
            est_mixmat  <<- read_ext('est_mixmat_log_log_scaled.rds')
            est_pcr <<- est_senesence <<- array(1, c(66, 2))
        }
    ),
    active = list()
)