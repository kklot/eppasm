#include "popClass.hpp"
#include "hivClass.hpp"
#include "artClass.hpp"

template<typename T>
class Model
{
public:
  StateSpace* const s;
  Parameters* const p;
  popC<T> pop;
  hivC<T> hivpop;
  artC<T> artpop;
public:
  Model(StateSpace* const ss, Parameters* const pp) :
    s(ss),
    p(pp),
    pop(s, p),
    hivpop(s, p),
    artpop(s, p)
    {}
public:
  void simulate(int t) {
    ++s->year;
    pop.update_views(t);
    hivpop.update_views(t);
    artpop.update_views(t);
    aging();
    death();
    migration();
    pop.update_fertile();
    if (s->MODEL != 0)  {
      infection_process();
      if (p->ic.eppmod == 2)
        pop.epp_disease_model_direct(hivpop, artpop);
    }
    adjust_pop();
    save_outputs();
  }

  void aging() 
  {
    pop.aging();
    pop.add_entrants();
    if (s->MODEL == 2)
      pop.sexual_debut();
    if (s->MODEL != 0) {
      pop.update_hiv_aging_prob();
      hivpop.aging(pop.hiv_aging_prob_);
      hivpop.add_entrants(pop.entrant_art_);
      if (s->MODEL == 2)
        hivpop.sexual_debut();
      if (s->year > s->tARTstart - 1) {
        artpop.aging(pop.hiv_aging_prob_);
        artpop.add_entrants(pop.entrant_art_);
        if (s->MODEL == 2)
          artpop.sexual_debut();
      }
    }
  }

  void death() {
    pop.deaths();
    if (s->MODEL != 0) {
      hivpop.deaths(pop.hiv_sx_prob);
      if (s->year > s->tARTstart - 1)
        artpop.deaths(pop.hiv_sx_prob);
    }
  }

  void migration() {
    pop.migration();
    if (s->MODEL != 0) {
      hivpop.migration(pop.hiv_mr_prob);
      if (s->year > s->tARTstart - 1)
        artpop.migration(pop.hiv_mr_prob);
    }
  }

  void adjust_pop() {
    if (p->dm.flag_popadjust) { // match target pop
      pop.adjust_pop();
      if (s->MODEL != 0) {
        hivpop.adjust_pop(pop.adj_prob);
        if (s->year >= s->tARTstart - 1)
          artpop.adjust_pop(pop.adj_prob);
      }
    }
  }

  void infection_process() 
  {
    for (int time_step = 0; time_step < s->steps_per_year; ++time_step) 
    {
      if (p->ic.eppmod != 2) { // != "directincid"
        pop.update_rvec(time_step);
        if (s->MIX)
          pop.infect_mix(hivpop, artpop, time_step);
        else
          pop.infect_spec(hivpop, artpop, time_step);
        pop.update_infection();
        hivpop.update_infection(pop.infections_);
      }
      if (p->ad.scale_cd4_mort && s->year >= s->tARTstart - 1)
        hivpop.scale_cd4_mort(artpop);
      hivpop.grad_progress(); // cd4 disease progression and mortality
      if (s->year >= s->tARTstart - 1)
        artpop.count_death();
      pop.remove_hiv_death(hivpop, artpop); // Remove hivdeaths from pop
      if (s->year >= s->tARTstart - 1) // ART initiation
        pop.epp_art_init(hivpop, artpop, time_step);
      hivpop.add_grad_to_pop();
    } // end time step
  }

  void save_outputs() {
    if (s->MODEL != 0) {
      pop.cal_prev_pregant(hivpop, artpop); // prevalence among pregnant women
      pop.save_prev_n_inc(); // save prevalence and incidence 15 to 49
    }
  }

};
