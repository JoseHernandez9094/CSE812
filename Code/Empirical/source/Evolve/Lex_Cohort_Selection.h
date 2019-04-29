#ifndef EMP_EVO_COHORT_SELECT_H
#define EMP_EVO_COHORT_SELECT_H

#include <map>
#include <functional>
#include <utility>

#include "../base/array.h"
#include "../base/assert.h"
#include "../base/vector.h"
#include "../base/macros.h"
#include "../meta/reflection.h"
#include "../tools/IndexMap.h"
#include "../tools/Random.h"
#include "../tools/vector_utils.h"

namespace emp {

  template<typename ORG> class World;

  /// ==SingleHostSelection== Selection picks a single individual to from a group of qualified applicants
  /// through a single fitness function
  /// world -> World we are evolving
  /// qulifieres -> id's of 
  /// fit_fun -> fitness function we are using to evaluate
  /// return -> a vector of ID's that have the highest fitness

  template<typename ORG>
  std::vector<size_t> SingleHostSelection(World<ORG> & world, const emp::vector<size_t> & qualifiers
                                                          , const std::function<double(ORG&)> & fit_fun)
  {
    std::map<double, std::vector<size_t>> scores;
    
    for(size_t id : qualifiers)
    {
      auto & org = world.GetOrg(id);
      double score = fit_fun(org);
      auto it = scores.find(score);

      //if found insert into existing vector
      if(it != scores.end())
      {
        scores[score].push_back(id);
      }
      //else add new vector to the map with score
      else
      {
        std::vector<size_t> temp;
        temp.push_back(id);
        scores.insert(std::pair<double, std::vector<size_t>>(score, temp));
      }
    }
    return scores.rbegin()->second;
  }

  /// ==CohortHostSelect== Selection picks a set of the most fit individuals from a certian
  /// cohort. It will use a lexicase selection method where we will keep looking 
  /// through all the fitness functions until we have a single winner. In the case
  /// that there is no single winner, we will randomly select one of the remaining
  /// at random
  /// world -> World we are evolving
  /// cohorts -> vector of cohorts
  /// fit_funs -> vector of fitness functions we are checking
  /// cohort_size -> size of an individual cohort
  /// cohort_total -> total number of cohorts

  template<typename ORG>
  void CohortHostSelect(World<ORG> & world, const emp::vector<emp::vector<size_t>> & cohorts, 
                                    emp::vector<std::function<double(ORG&)>> fit_funs,
                                    const size_t & cohort_size, const size_t & cohort_total, emp::Ptr<emp::Random> rng) 
  {
    std::vector<size_t> winners;
    
    for(std::vector<size_t> cohort: cohorts)
    {
      size_t win_cnt = 0;
      
      while(win_cnt != cohort_size)
      {
        size_t fun_ptr = 0;
        emp::Shuffle(*rng, fit_funs);
        std::vector<size_t> qualifiers = SingleHostSelection(world, cohort, fit_funs[fun_ptr]);
        ++fun_ptr;

        while(qualifiers.size() > 1 && fun_ptr != fit_funs.size())
        {
          qualifiers = SingleHostSelection(world, qualifiers, fit_funs[fun_ptr]);
          ++fun_ptr;
        }

        //If number of qualified ID's is just 1
        if(qualifiers.size() == 1)
        {
          winners.push_back(qualifiers[0]);
        }
        //Qualified ID's > 1 and fun_ptr == fit_funs.size
        else
        {
          std::vector<size_t> choice = emp::Choose(*rng, qualifiers.size(), 1);
          winners.push_back(qualifiers[choice[0]]);
        }
        
        ++win_cnt;
      }
    }

    if(winners.size() != (cohort_size * cohort_total))
    {
      std::cout << "Cohort_Selection winners vector not correct size!" << std::endl;
      exit(-1);
    }

    for(size_t id : winners)
    {
      world.DoBirth(world.GetGenomeAt(id), id, 1);
    }
  }
}

#endif
