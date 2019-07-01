#include "Classes.hpp"

void artC::aging(const boost2D& ag_prob, Views& v, const StateSpace& s) {
  std::memcpy(at_this, at_prev, N*sizeof(double));
    if (s.MODEL == 2)
      std::memcpy(at_this_db, at_prev_db, N*sizeof(double));
  double nARTup;
  int agr_size = s.hDS * s.hTS;
  int k = 0, offset = N / s.NG,  offset_k = s.hAG;
  for (int i = 0; i < agr_size * (s.hAG-1); i += agr_size) {
    for (int j = 0; j < agr_size; j++) {
      nARTup = *(at_prev + i+j) * *(ag_prob.data() + k);
      *(at_this + i+j)          -= nARTup;
      *(at_this + i+j+agr_size) += nARTup;
      nARTup = *(at_prev + offset+i+j) * *(ag_prob.data() + k + offset_k);
      *(at_this + offset+i+j)          -= nARTup;
      *(at_this + offset+i+j+agr_size) += nARTup;
    }
    ++k;
  }

  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG-1; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          if (s.MODEL == 2 && agr < s.hDB - 1) {
            nARTup = data_db[s.year-1][sex][agr][cd4][dur] * ag_prob[sex][agr];
            data_db[s.year][sex][agr][cd4][dur]   -= nARTup;
            data_db[s.year][sex][agr+1][cd4][dur] += nARTup;
          }
        }
}

void artC::add_entrants(const dvec& artYesNo, Views& v, const Parameters& p,
                        const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int cd4 = 0; cd4 < s.hDS; cd4++)
      for (int dur = 0; dur < s.hTS; dur++) {
        double add =
          p.ph.paedsurv_artcd4dist[s.year][sex][cd4][dur] * artYesNo[sex];
        if (s.MODEL == 1)
          v.now_art[sex][0][cd4][dur] += add;
        if (s.MODEL == 2) // add to virgin then debut
          data_db[s.year][sex][0][cd4][dur] += add;
      }
}

void artC::sexual_debut(Views& v, const Parameters& p, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hDB; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          double n_db = data_db[s.year][sex][agr][cd4][dur] * p.ic.db_pr[sex][agr];
          v.now_art[sex][agr][cd4][dur]       += n_db;
          data_db[s.year][sex][agr][cd4][dur] -= n_db;
        }
}

void artC::deaths(const boost2D& survival_pr, Views& v, const StateSpace& s) {
  const double * at_sx = survival_pr.data(); int agr = 0;
  for (int i = 0; i < N; i += s.hTS * s.hDS) {
    for (int j = 0; j < s.hTS * s.hDS; j++) {
      *(at_this + i+j) *= *(at_sx + agr);
      if (s.MODEL == 2)
        *(at_this_db + i+j) *= *(at_sx + agr);
    }
    ++agr;
  }
}

void artC::migration(const boost2D& migration_pr, Views& v, const StateSpace& s) {
  const double * at_mr = migration_pr.data(); int agr = 0;
  for (int i = 0; i < N; i += s.hTS * s.hDS) {
    for (int j = 0; j < s.hTS * s.hDS; j++) {
      *(at_this + i+j) *= *(at_mr + agr);
      if (s.MODEL == 2)
        *(at_this_db + i+j) *= *(at_mr + agr);
    }
    ++agr;
  }
}

void artC::grad_progress(Views& v, const StateSpace& s) {
  for (int i = 0; i < N; i++) {
    *(gradART.data() + i) = 0; // reset gradient
    *(gradART.data() + i) -= *(death_.data() + i);
    if (s.MODEL == 2) {
      *(gradART_db.data() + i) = 0; // reset gradient
      *(gradART_db.data() + i) -= *(death_db_.data() + i);
    }
  }
  double art_up;
  for (int i = 0; i < N; i += s.hTS) {
    for (int j = 0; j < s.hTS - 1; j++) {
      art_up = 2. * *(at_this +i+j);
      *(gradART.data() + i+j)   -= art_up;
      *(gradART.data() + i+j+1) += art_up;
      if (s.MODEL == 2) {
        art_up = 2. * *(at_this_db +i+j);
        *(gradART_db.data() + i+j)   -= art_up;
        *(gradART_db.data() + i+j+1) += art_up;
      }
    }
  }
}

void artC::art_dropout(hivC& hivpop, Views& v,
                       const Parameters& p,
                       const StateSpace& s) {
  double n_dropout, p_dropout = p.ad.art_dropout[s.year]; int k = 0;
  for (int i = 0; i < N; i += s.hTS) {
    for (int j = 0; j < s.hTS; j++) {
      n_dropout = *(at_this + i+j)  * p_dropout;
      *(hivpop.grad.data()  +   k) += n_dropout;
      *(gradART.data() + i+j)      -= n_dropout;
      if (s.MODEL == 2) {
        n_dropout = *(at_this_db + i+j) * p_dropout;
        *(hivpop.grad_db.data()  +   k) += n_dropout;
        *(gradART_db.data()      + i+j) -= n_dropout;        
      }
    }
    ++k;
  }
}

void artC::update_current_on_art(Views& v, const StateSpace& s) {
  art_by_sex_.assign(s.NG, .0); int sex = 0; // reset when call
  for (int i = 0; i < N; i += N/2) {
    for (int j = 0; j < N/2; j++) {
      art_by_sex_[sex] += *(at_this +i+j) + *(gradART.data() + i+j) * s.DT;
      if (s.MODEL == 2)
        art_by_sex_[sex] += *(at_this_db +i+j) + *(gradART_db.data() + i+j) * s.DT;
    }
    ++sex;
  }
}

void artC::grad_init(const boost3D& artinit, Views& v, const StateSpace& s) {
  int j = 0;
  for (int i = 0; i < N; i += s.hTS) {
    gradART.data()[i] += artinit.data()[j] / s.DT;
    ++j;
  }
  for (int i = 0; i < N; i++)
    *(at_this + i) += *(gradART.data() + i) * s.DT;
}

void artC::grad_db_init(const boost3D& artinit_db, const StateSpace& s) {
  int j = 0;
  for (int i = 0; i < N; i += s.hTS) {
    gradART_db.data()[i] += artinit_db.data()[j] / s.DT;
    ++j;
  }
  for (int i = 0; i < N; i++)
    *(at_this_db + i) += *(gradART_db.data() + i) * s.DT;
}

void artC::adjust_pop(const boost2D& adj_prob, Views& v, const StateSpace& s) {
  int agr = 0;
  for (int i = 0; i < N; i += s.hTS * s.hDS) {
    for (int j = 0; j < s.hTS * s.hDS; j++) {
      *(at_this + i+j) *= *(adj_prob.data() + agr);
      if (s.MODEL == 2)
        *(at_this_db + i+j) *= *(adj_prob.data() + agr);
    }
    ++agr;
  }
}

void artC::count_death(Views& v, const Parameters& p, const StateSpace& s) {
  double * at_death = death_.data(), * at_death_db = death_db_.data();
  const double * at_amx = p.nh.art_mort.data(),
               * at_mxrr = p.nh.artmx_timerr.data();
  for (int i = 0; i < N; i += s.hTS)
    for (int j = 0; j < s.hTS; j++) {
      *(at_death + i+j) = *(at_this + i+j) * *(at_amx + i+j) * *(at_mxrr + j);
      if (s.MODEL == 2)
        *(at_death_db + i+j) = *(at_this_db + i+j) * *(at_amx + i+j) * *(at_mxrr + j);
    }
}