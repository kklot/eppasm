// forward declare
template<typename T> epp::matrix<T> getVmat(const Parameters& p, int sex);

template<typename T>
epp::VectorX<T> NGM(const Parameters& p, const StateSpace& s) 
{
  int N  = p.dm.basepop.dimension(0);
  epp::matrix<T> F_form(8,8), F_scale(8,8);
  F_form.setZero();
  F_scale.setZero();
  F_form.chip(0, 0) = F_form.chip(0, 0) + 1.;
  F_scale(0,0) = 1.;

  epp::matrix<T> scale_mat(N, N);
  scale_mat.setConstant(p.ic.stage0_kappa);
  epp::matrix<T> scales = epp::kronecker<T>(F_scale, scale_mat);
  scales = (scales == 0.).select(scales.constant(1), scales);

  epp::vector<T> 
    _m = p.ic.rvec(0) * p.ic.est_pcr.chip(s.M, 1) * p.dm.basepop.chip(s.M, 1),
    _f = p.ic.rvec(0) * p.ic.est_pcr.chip(s.F, 1) * p.dm.basepop.chip(s.F, 1);

  epp::matrix<T> 
    // FOI
    Ffm0i_ = epp::sweep<T>(p.ic.mixmat.chip(s.M, 2), 0, _m),
    Fmf0i_ = epp::sweep<T>(p.ic.mixmat.chip(s.F, 2), 0, _f),
    Ffm0i  = epp::sweep<T>(Ffm0i_, 1, p.dm.basepop.chip(s.F, 1), "/"),
    Fmf0i  = epp::sweep<T>(Fmf0i_, 1, p.dm.basepop.chip(s.M, 1), "/"),
    Ffm    = epp::kronecker<T>(F_form, Ffm0i)  * scales,
    Fmf    = epp::kronecker<T>(F_form, Fmf0i)  * scales,
    Vm     = getVmat<T>(p, s.M),
    Vf     = getVmat<T>(p, s.F),
    FF     = epp::kronecker<T>(epp::matrix_single_entry<T>(4,1,0), Fmf) + 
             epp::kronecker<T>(epp::matrix_single_entry<T>(4,0,1), Ffm),
    VV     = epp::kronecker<T>(epp::matrix_single_entry<T>(4,0,0), Vm) + 
             epp::kronecker<T>(epp::matrix_single_entry<T>(4,1,1), Vf),
    Jac    = FF - VV;
  epp::map_of_MatrixX<T> Jac_map(&Jac(0), Jac.dimension(0), Jac.dimension(1));
  epp::VectorX<T> leading_ev = epp::power_method<T>(Jac_map);
  leading_ev /= leading_ev.sum();
  return leading_ev;
}

template<typename T>
epp::matrix<T> getVmat(const Parameters& p, int sex)
{
  int N = p.dm.basepop.dimension(0);
  T DT = 1., alpha = 1. - exp(-1./(p.ic.stage0_time * 12. * DT));
  epp::matrix<T>
    rhos = p.nh.cd4_prog_xp.chip(sex, 2),
    ini = p.nh.cd4_initdist_xp.chip(sex, 2);
  epp::vector<T> death = 1. - p.dm.Sx.chip(0, 2).chip(sex, 1);

  epp::matrix<T> Vjj = epp::Diag<T>(1. + death); // core diag matrix
  for (int i = 1; i < N; i++) // its lower diag
    Vjj(i,i-1) = - death(i-1);

  epp::matrix<T>
    // ageing and death
    V00 = Vjj + epp::Diag<T>(alpha, N),
    V11 = Vjj + epp::Diag<T>(rhos.chip(1-1, 1)),
    V22 = Vjj + epp::Diag<T>(rhos.chip(2-1, 1)),
    V33 = Vjj + epp::Diag<T>(rhos.chip(3-1, 1)),
    V44 = Vjj + epp::Diag<T>(rhos.chip(4-1, 1)),
    V55 = Vjj + epp::Diag<T>(rhos.chip(5-1, 1)),
    V66 = Vjj + epp::Diag<T>(rhos.chip(6-1, 1)), 
    V77 = Vjj,
    // intial infection to other stages
    V10 = epp::Diag<T>(-alpha * ini.chip(1-1, 1)),
    V20 = epp::Diag<T>(-alpha * ini.chip(2-1, 1)),
    V30 = epp::Diag<T>(-alpha * ini.chip(3-1, 1)),
    V40 = epp::Diag<T>(-alpha * ini.chip(4-1, 1)),
    V50 = epp::Diag<T>(-alpha * ini.chip(5-1, 1)),
    V60 = epp::Diag<T>(-alpha * ini.chip(6-1, 1)),
    V70 = epp::Diag<T>(-alpha * ini.chip(7-1, 1)), 
    // stages progression
    V21 = epp::Diag<T>(rhos.chip(1-1, 1)),
    V32 = epp::Diag<T>(rhos.chip(2-1, 1)),
    V43 = epp::Diag<T>(rhos.chip(3-1, 1)),
    V54 = epp::Diag<T>(rhos.chip(4-1, 1)),
    V65 = epp::Diag<T>(rhos.chip(5-1, 1)),
    V76 = epp::Diag<T>(rhos.chip(6-1, 1)),
    // the V matrix of one sex
    aV =
      // diag
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 0, 0), V00) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 1, 1), V11) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 2, 2), V22) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 3, 3), V33) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 4, 4), V44) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 5, 5), V55) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 6, 6), V66) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 7, 7), V77) +
      // left most
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 1, 0), V10) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 2, 0), V20) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 3, 0), V30) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 4, 0), V40) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 5, 0), V50) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 6, 0), V60) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 7, 0), V70) +
      // lower diag
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 2, 1), V21) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 3, 2), V32) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 4, 3), V43) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 5, 4), V54) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 6, 5), V65) +
      epp::kronecker<T>(epp::matrix_single_entry<T>(8, 7, 6), V76);
  return aV;
}
