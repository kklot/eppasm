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
