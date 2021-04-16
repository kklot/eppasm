#' @useDynLib eppasm eppasmC eppasmOOpp
foo <- function() {}

#' Using power method to obtain the dominance eigenvector from infection matrix
#' 
#' For large (sparse) matrix this speeds up up to 100x compare to compute all 
#' eigenvectors using base R's eigen.
#' 
#' @param x the Jacobian of the infection subsystem, or F+V in Next Generation
#' Matrix, will be automatically converted to sparse
#' @param error error tolerance, using 1e-9 would make most result comparable to base R eigen's result in most cases with a slightly increase in time.
#' @param R0 initial guess of the dominance eigenvalue
#' @param max_iter if the algorithm does not converged to the error tolerance, 
#' do a hard stop at this iteration (with a warning).
#' @useDynLib eppasm power_method_symbol
#' @export
power_method <- function(x, error = 1e-6, R0 = 1, max_iter = 1000L) {
    .Call("power_method_symbol", 
          as.matrix(x), nrow(x), ncol(x), error, R0, max_iter)
}
