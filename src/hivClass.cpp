#include "Classes.hpp"

void hivC::aging(const boost2D& ag_prob, Views& v, const StateSpace& s) {
  std::memcpy(at_this, at_prev, N * sizeof(double));
    if (s.MODEL == 2)
      std::memcpy(at_this_db, at_prev_db, N * sizeof(double));
  double nHup;
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG - 1; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        nHup = v.pre_hiv[sex][agr][cd4] * ag_prob[sex][agr];
        v.now_hiv[sex][agr][cd4]   -= nHup;
        v.now_hiv[sex][agr+1][cd4] += nHup;
        if (s.MODEL == 2 && agr < s.hDB - 1) {
          nHup = data_db[s.year-1][sex][agr][cd4] * ag_prob[sex][agr];
          data_db[s.year][sex][agr][cd4]   -= nHup;
          data_db[s.year][sex][agr+1][cd4] += nHup;
        }
      }
}

void hivC::add_entrants(const dvec& artYesNo, Views& v,
                        const Parameters& p, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int cd4 = 0; cd4 < s.hDS; cd4++) {
      double add = p.ph.paedsurv_cd4dist[s.year][sex][cd4] * artYesNo[sex+2];
      if (s.MODEL == 1)
        v.now_hiv[sex][0][cd4] += add;
      if (s.MODEL == 2) // add to virgin then debut
        data_db[s.year][sex][0][cd4] += add;
    }
}

void hivC::sexual_debut(Views& v, const Parameters& p, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int adb = 0; adb < s.hDB; adb++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        double n_db = data_db[s.year][sex][adb][cd4] * p.ic.db_pr[sex][adb];
        v.now_hiv[sex][adb][cd4]       += n_db;
        data_db[s.year][sex][adb][cd4] -= n_db;
      }
}

void hivC::deaths (const boost2D& survival_pr, Views& v, const StateSpace& s) {
  const double * at_surv = survival_pr.data(); int agr = 0;
  for (int i = 0; i < N; i += s.hDS) {
    for (int j = 0; j < s.hDS; j++) {
      *(at_this + i+j) *= *(at_surv + agr);
      if (s.MODEL == 2)
        *(at_this_db + i+j) *= *(at_surv + agr);
    } 
    ++agr;
  }
}

void hivC::migration (const boost2D& migration_pr, Views& v, const StateSpace& s) {
  const double * at_mr = migration_pr.data(); int agr = 0;
  for (int i = 0; i < N; i += s.hDS) {
    for (int j = 0; j < s.hDS; j++) {
      *(at_this + i+j) *= *(at_mr + agr);
      if (s.MODEL == 2)
        *(at_this_db + i+j) *= *(at_mr + agr);
    } 
    ++agr;
  }
}

void hivC::update_infection (const boost2D& new_infect,
                             const Parameters& p, const StateSpace& s) {
  int j = 0, span = 0;
  for (int i = 0; i < s.NG * s.pAG; i += span) {
    infect_by_agrp_.data()[j] = 0;
    span = (j >= s.hAG) ? s.h_ag_span[j - s.hAG] : s.h_ag_span[j];
    for (int k = 0; k < span; ++k)
      infect_by_agrp_.data()[j] += new_infect.data()[i+k];
    ++j;
  }
  const double * at_cd4_initdist = p.nh.cd4_initdist.data(); int agr = 0;
  for (int i = 0; i < N; i += s.hDS) {
    for (int j = 0; j < s.hDS; j++) {
      grad.data()[i+j]  = 0.; // reset every time step
      grad.data()[i+j] += infect_by_agrp_.data()[agr] * at_cd4_initdist[i+j];
    }
    ++agr;
  }
}

void hivC::scale_cd4_mort(artC& artpop, Views& v,
                           const Parameters& p, const StateSpace& s) {
  double num, den = 0;
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        num = v.now_hiv[sex][agr][cd4] + data_db[s.year][sex][agr][cd4];
        for (int dur = 0; dur < s.hTS; dur++)
          den += v.now_art[sex][agr][cd4][dur] +
                 artpop.data_db[s.year][sex][agr][cd4][dur];
        num = (num + den == 0.0) ? 1 : num / (num + den);
        cd4_mort_[sex][agr][cd4] = num * p.nh.cd4_mort[sex][agr][cd4];
        den = 0;
      }
}

void hivC::grad_progress(Views& v, const Parameters& p, const StateSpace& s) {
  // HIV death
  for (int i = 0; i < N; i++) {
    *(death_.data() + i)  = *(at_this + i) * *(cd4_mort_.data() + i);
    if (p.ic.eppmod == 2) 
      *(grad.data() + i)  = 0.; // reset in direct incidence input model
    *(grad.data()   + i) -= *(death_.data() + i);
    if (s.MODEL == 2) {
      *(death_db_.data() + i)  = *(at_this_db + i) * *(cd4_mort_.data() + i);
      *(grad_db.data()   + i)  = 0.; // reset, 1st time grad_db is used
      *(grad_db.data()   + i) -= *(death_db_.data() + i);
    }
  }
  
  // remove cd4 stage progression (untreated)
  int k = 0; const double * at_cd4_prog = p.nh.cd4_prog.data();
  for (int i = 0; i < N; i += s.hDS) {
    for (int j = 0; j < s.hDS - 1; j++) {
      double nHup = *(at_cd4_prog + k+j) * *(at_this + i+j);
      *(grad.data() + i+j)     -= nHup;
      *(grad.data() + i+j+1)   += nHup;
      if (s.MODEL == 2) {
        nHup = *(at_cd4_prog + k+j) * *(at_this_db + i+j);
        *(grad_db.data() + i+j)     -= nHup;
        *(grad_db.data() + i+j+1)   += nHup;
      }
    }
    k += s.hDS - 1; // cd4_prog is (hDS-1)xhAGxNG not hDSxhAGxNG
  }
}

void hivC::distribute_artinit (boost3D& artinit, artC& artpop,
                               Views& v, const StateSpace& s) {
    double debut_now, all_hivpop, pr_weight_db, n_artinit_db;
    boost3D artinit_db(extents[s.NG][s.hAG][s.hDS]);
    for (int sex = 0; sex < s.NG; sex++)
      for (int agr = 0; agr < s.hAG; agr++)
        for (int cd4 = 0; cd4 < s.hDS; cd4++) {
          debut_now =
            data_db[s.year][sex][agr][cd4] + s.DT * grad_db[sex][agr][cd4];
          all_hivpop =
            (v.now_hiv[sex][agr][cd4] + s.DT * grad[sex][agr][cd4]) + debut_now;
          if (artinit[sex][agr][cd4] > all_hivpop)
            artinit[sex][agr][cd4]   = all_hivpop;
          pr_weight_db               = debut_now / all_hivpop;
          n_artinit_db               = artinit[sex][agr][cd4] * pr_weight_db;
          artinit_db[sex][agr][cd4]  = n_artinit_db;
          artinit[sex][agr][cd4]    -= n_artinit_db;
          grad_db[sex][agr][cd4]    -= n_artinit_db / s.DT;
        }
    artpop.grad_db_init(artinit_db, s);
}

void hivC::add_grad_to_pop (Views& v, const StateSpace& s) {
  for (int i = 0; i < N; ++i) {
    *(at_this + i) += s.DT * *(grad.data() + i);
    if (s.MODEL == 2)
      *(at_this_db + i) += s.DT * *(grad_db.data() + i);
  }
}

void hivC::adjust_pop(const boost2D& adj_prob, Views& v, const StateSpace& s) {
  int agr = 0;
  for (int i = 0; i < N; i += s.hDS) {
    for (int j = 0; j < s.hDS; j++) {
      *(at_this + i+j) *= *(adj_prob.data() + agr);
      if (s.MODEL == 2)
        *(at_this_db + i+j) *= *(adj_prob.data() + agr);
    } 
    ++agr;
  }
}