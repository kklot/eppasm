#' Obtain "V" matrix for the NGM
#' 
#' @param p model parameters (fp)
#' @param sex male=1, female=2
getVmat = function(p, sex) {
  lower_diag = function(x) cbind(2:nrow(x), 1:(nrow(x)-1))
  N = nrow(p$basepop)
  death = 1 - p$Sx[,sex,1]
  Vjj = Matrix::Diagonal(x=death + 1) 
  Vjj[ lower_diag(Vjj) ] = -death[1:(N-1)]

  DT = 1                                # doing 1 year time-step
  alpha = 1 - exp(-1/(p$stage0_time * 12 * DT))
  rhos = p$cd4_prog[,,sex] %>% apply(1, rep, times=p$ss$h.ag.span)
  
  V00 = Vjj + Matrix::Diagonal(N, alpha)

  V11 = Vjj + Matrix::Diagonal(x = rhos[, 1])
  V22 = Vjj + Matrix::Diagonal(x = rhos[, 2])
  V33 = Vjj + Matrix::Diagonal(x = rhos[, 3])
  V44 = Vjj + Matrix::Diagonal(x = rhos[, 4])
  V55 = Vjj + Matrix::Diagonal(x = rhos[, 5])
  V66 = Vjj + Matrix::Diagonal(x = rhos[, 6])
  
  V77 = Vjj

  ini = p$cd4_initdist[,,sex] %>% apply(1, rep, times=p$ss$h.ag.span)

  V10 = Matrix::Diagonal(x = -alpha * ini[, 1])
  V20 = Matrix::Diagonal(x = -alpha * ini[, 2])
  V30 = Matrix::Diagonal(x = -alpha * ini[, 3])
  V40 = Matrix::Diagonal(x = -alpha * ini[, 4])
  V50 = Matrix::Diagonal(x = -alpha * ini[, 5])
  V60 = Matrix::Diagonal(x = -alpha * ini[, 6])
  V70 = Matrix::Diagonal(x = -alpha * ini[, 7])

  V21 = Matrix::Diagonal(x=rhos[, 1])
  V32 = Matrix::Diagonal(x=rhos[, 2])
  V43 = Matrix::Diagonal(x=rhos[, 3])
  V54 = Matrix::Diagonal(x=rhos[, 4])
  V65 = Matrix::Diagonal(x=rhos[, 5])
  V76 = Matrix::Diagonal(x=rhos[, 6])

  diags = list(V00, V11, V22, V33, V44, V55, V66, V77)

  Matrix::.bdiag(diags) +
    kronecker(Matrix::Matrix(c(rep(0,1), 1, rep(0,8*8-1-1)), 8), V10) +
    kronecker(Matrix::Matrix(c(rep(0,2), 1, rep(0,8*8-1-2)), 8), V20) +
    kronecker(Matrix::Matrix(c(rep(0,3), 1, rep(0,8*8-1-3)), 8), V30) +
    kronecker(Matrix::Matrix(c(rep(0,4), 1, rep(0,8*8-1-4)), 8), V40) +
    kronecker(Matrix::Matrix(c(rep(0,5), 1, rep(0,8*8-1-5)), 8), V50) +
    kronecker(Matrix::Matrix(c(rep(0,6), 1, rep(0,8*8-1-6)), 8), V60) +
    kronecker(Matrix::Matrix(c(rep(0,7), 1, rep(0,8*8-1-7)), 8), V70) +
    kronecker(Matrix::Matrix(c(rep(0,8*1 + 2 + 0), 1, rep(0,8*8 - 1 - 8*1 - 2 - 0)), 8), V21) +
    kronecker(Matrix::Matrix(c(rep(0,8*2 + 2 + 1), 1, rep(0,8*8 - 1 - 8*2 - 2 - 1)), 8), V32) +
    kronecker(Matrix::Matrix(c(rep(0,8*3 + 2 + 2), 1, rep(0,8*8 - 1 - 8*3 - 2 - 2)), 8), V43) +
    kronecker(Matrix::Matrix(c(rep(0,8*4 + 2 + 3), 1, rep(0,8*8 - 1 - 8*4 - 2 - 3)), 8), V54) +
    kronecker(Matrix::Matrix(c(rep(0,8*5 + 2 + 4), 1, rep(0,8*8 - 1 - 8*5 - 2 - 4)), 8), V65) +
    kronecker(Matrix::Matrix(c(rep(0,8*6 + 2 + 5), 1, rep(0,8*8 - 1 - 8*6 - 2 - 5)), 8), V76) 
}

#' Next generation matrix for the sexual mixing model
#'
#' this is simplified as the absolute size is not important to the relative size
#'
#' @param p list of sexual behaviors model's parameters
#' @param kappa increase of infectiousness in stage 0 vs 1-7
#' @export 
domimance_vector = function(p, stage0_kappa = 10, scale = TRUE, full = FALSE, version = 'C'){
  p$stage0_kappa = stage0_kappa
	if (version == 'C') 
  {
    V0 = .Call("domimance_vector_symbol", prepare_fp_for_Cpp(p)) #  scaled in C
  }
  else {
    N = nrow(p$basepop)
    # F matrix
    F_form = F_scale = Matrix::Matrix(0, 8, 8)
    F_form[1, ] = 1
    F_scale[1,1] = 1
  
    scale_mat = Matrix::Matrix(stage0_kappa, N, N)
    scales = kronecker(F_scale, scale_mat)
    scales[scales==0] = 1
  
    Ffm0i = (p$rvec[1] * p$est_pcr[,1] * sweep(p$mixmat[,,1], 1, p$basepop[, 1], "*")) %>%
      sweep(2, p$basepop[, 2], "/") %>% as("sparseMatrix")
    Fmf0i = (p$rvec[1] * p$est_pcr[,2] * sweep(p$mixmat[,,2], 1, p$basepop[, 2], "*")) %>%
      sweep(2, p$basepop[, 1], "/") %>% as("sparseMatrix")
  
    Ffm = kronecker(F_form, Ffm0i) * scales
    Fmf = kronecker(F_form, Fmf0i) * scales
    FF = kronecker(Matrix::Matrix(c(0,1,0,0), 2), Fmf) + kronecker(Matrix::Matrix(c(0,0,1,0), 2), Ffm)
  
    # V matrix
    Vm = getVmat(p, 1)
    Vf = getVmat(p, 2)
    VV = kronecker(Matrix::Matrix(c(1,0,0,0), 2), Vm, sparse=T) + 
         kronecker(Matrix::Matrix(c(0,0,0,1), 2), Vf, sparse=T)
    Jac = FF - VV
    if (full) {
      ei = eigen(Jac)
      return(list(Jac=Jac, R0 = abs(Re(ei$values[1])), V0 = abs(Re(ei$vectors[,1]))))
    }  
    V0  = power_method(Jac)
    if (scale)
      V0  = V0 / sum(V0)
  }
  # reshape to model dimensions
  V0  <- aperm(array(V0, c(p$ss$pAG, p$ss$hDS + 1, p$ss$NG)), c(1,3,2))
  # with approxmation there are artifacts in the last ages from power_method, to
  # remove it the error need to reduce to somewhere like 1e-12 and cost time 
  # while the general shape of the vector is not different, so set to zero here
  V0[a2i(70:80),,] = 0 # age x sex x (ds + 1)
  V0

}
