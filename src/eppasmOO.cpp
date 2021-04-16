#include "Rinterface.hpp"

extern "C" SEXP eppasmOOpp(SEXP fp) {
  StateSpace s(fp); // read only state-space
  Parameters p(fp); // read only parameters
  Model<double> model(&s, &p); 
  for (int i = 1; i < s.SIM_YEARS; ++i)
    model.simulate(i);
  outputSEXP<double> O(s, &model);  // create all SEXP for output
  return O.pop;
}

extern "C" SEXP power_method_symbol(SEXP Jac, SEXP nrow, SEXP ncol, 
                                    SEXP e, SEXP R0, SEXP max_iter) {
  Eigen::Map<Eigen::MatrixXd> aJac(REAL(Jac), *INTEGER(nrow), *INTEGER(ncol));
  Eigen::VectorXd leading_ev = 
    epp::power_method<double>(aJac, *REAL(e), *REAL(R0), *INTEGER(max_iter));
  SEXP r_ev = PROTECT(NEW_NUMERIC(*INTEGER(nrow)));
  memcpy(REAL(r_ev), leading_ev.data(), sizeof(double) * *INTEGER(nrow));
  UNPROTECT(1);
  return r_ev;
}
