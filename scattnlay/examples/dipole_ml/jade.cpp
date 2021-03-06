///
/// @file   jade.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Thu Aug 15 19:26:53 2013
/// @copyright 2013 Ladutenko Konstantin
/// @section LICENSE
/// This file is part of JADE++.
///
/// JADE++ is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// JADE++ is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with JADE++.  If not, see <http://www.gnu.org/licenses/>.
///
/// @brief JADE++ is a free (GPLv3+) high performance implementation of
/// adaptive differential evolution optimization algorithm from
/// Jingqiao Zhang and Arthur C. Sanderson book 'Adaptive Differential
/// Evolution. A Robust Approach to Multimodal Problem Optimization'
/// Springer, 2009.  Crossover rate was patched according to PMCRADE
/// approach supposed by Jie Li, Wujie Zhu, Mengjun Zhou, and Hua Wang
/// in 'Power Mean Based Crossover Rate Adaptive Differential
/// Evolution' in H. Deng et al. (Eds.): AICI 2011, Part II, LNAI
/// 7003, pp. 34–41, 2011
#include "./jade.h"
#include <mpi.h>
#include <random>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <map>
#include <iterator>
#include <string>
#include <vector>

namespace jade {
  /// @todo Replace all simple kError returns with something meangfull. TODO change kError to throw exception
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::Selection(std::vector<double> crossover_u, long i)  {
    bool is_evaluated = false;
    double f_current = 0.0;
    for (auto f : evaluated_fitness_for_current_vectors_) {
      if (f.second == i) {
        f_current = f.first;
        is_evaluated = true;
      }
    }  // end of searching of pre-evaluated fitness for current individual
    if (!is_evaluated) error_status_ = kError;
    double f_best = evaluated_fitness_for_current_vectors_.front().first;
    double f_crossover_u = FitnessFunction(crossover_u);
    bool is_success = f_crossover_u > f_current
        || f_crossover_u == f_best;  //Selected for maxima search
    if (is_find_minimum_) is_success = !is_success;
    if (!is_success) {
      // Case of current x and f were new for current generation.
      x_vectors_next_generation_[i] = x_vectors_current_[i];
      for (auto &f : evaluated_fitness_for_next_generation_) {
        if (f.second == i) f.first = f_current;
      }  // end of saving old fitness value in new generation.
    } else {  // if is_success == true
      x_vectors_next_generation_[i] = crossover_u;
      for (auto &f : evaluated_fitness_for_next_generation_) 
        if (f.second == i) f.first = f_crossover_u;
      to_be_archived_best_A_.push_back(x_vectors_current_[i]);
      successful_mutation_parameters_S_F_.push_back(mutation_F_[i]);
      successful_crossover_parameters_S_CR_.push_back(crossover_CR_[i]);
      // if (process_rank_ == kOutput)
      //   printf("n%li f_new=%4.2f\n",i,f_crossover_u);
      //PrintSingleVector(crossover_u);
    }  // end of dealing with success crossover
    return kDone;
  } // end of int SubPopulation::Selection(std::vector<double> crossover_u);
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::ArchiveCleanUp() {
    auto archived_last=archived_best_A_.end();
    archived_best_A_.splice(archived_last, to_be_archived_best_A_);
    auto size_A = archived_best_A_.size();
    long initial_diff = size_A - subpopulation_;
    // if (process_rank_ == kOutput)
    //   printf("diff = %li size_A=%li subpop=%li \n ", initial_diff, size_A, subpopulation_);
    // if (initial_diff < 1) return kDone;
    // if (process_rank_ == kOutput)
    //   printf("diff = %li size_A=%li \n ", initial_diff, size_A);
    for (long i = 0; i < initial_diff; ++i) {
      long index_to_remove = randint(0, size_A - 1);
      auto element_A = archived_best_A_.begin();
      std::advance(element_A, index_to_remove);
      archived_best_A_.erase(element_A);
      --size_A;
    }
    const unsigned long new_size_A = archived_best_A_.size();
    if (new_size_A > subpopulation_) error_status_ = kError; //TODO change kError to throw exception
    return kDone;
  } // end of int SubPopulation:: ArchiveCleanUp();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::Adaption()  {    
    long elements = 0;
    double sum = 0.0;
    for (auto CR : successful_crossover_parameters_S_CR_) {
      sum += CR;
      ++elements;
    }  // end of collecting data for mean CR
    if (elements == 0) return kDone;
    // Should never be reached for symmetric filling of S_CR and S_F.
    if (successful_mutation_parameters_S_F_.size() == 0) return kDone;  
    double mean_a_CR = sum / static_cast<double>(elements);
    // Original JADE adaption of mu_CR.
    if (!isPMCRADE_) {
      adaptor_crossover_mu_CR_ =
        (1 - adaptation_frequency_c_) *  adaptor_crossover_mu_CR_
        + adaptation_frequency_c_ * mean_a_CR;
    } else {
      // PMCRADE patch for mu_CR
      double std_S_CR = 0;
      for (auto CR : successful_crossover_parameters_S_CR_)
        std_S_CR += pow2(CR - mean_a_CR);
      std_S_CR = sqrt(std_S_CR / static_cast<double>(elements));
      const double PMCRADE_const = 0.07;
      if (std_S_CR < PMCRADE_const) {
        adaptor_crossover_mu_CR_ =
          (1 - adaptation_frequency_c_) *  adaptor_crossover_mu_CR_
          + adaptation_frequency_c_ * mean_a_CR;
      } else {
        double mean_pow2_CR = 0.0;
        for (auto CR : successful_crossover_parameters_S_CR_)
          mean_pow2_CR += pow2(CR);
        mean_pow2_CR = sqrt(mean_pow2_CR/static_cast<double>(elements));
        adaptor_crossover_mu_CR_ =
          (1 - adaptation_frequency_c_) *  adaptor_crossover_mu_CR_
          + adaptation_frequency_c_ * mean_pow2_CR;
      }
      // end of PMCRADE patch
    }  // end of if isPMCRADE_
    double sum_F = 0.0, sum_F2 = 0.0;
    for (auto F : successful_mutation_parameters_S_F_) {
      sum_F += F;
      sum_F2 += F*F;
    }  // end of collection data for Lehmer mean F
    double mean_l_F = sum_F2 / sum_F;
    adaptor_mutation_mu_F_ =
      (1 - adaptation_frequency_c_) *  adaptor_mutation_mu_F_
      + adaptation_frequency_c_ * mean_l_F;
    return kDone;
  } // end of int SubPopulation::Adaption();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::RunOptimization() {
    //if (process_rank_ == kOutput) printf("Start optimization..\n");
    if (error_status_) return error_status_;
    adaptor_mutation_mu_F_ = 0.5;
    adaptor_crossover_mu_CR_ = 0.5;
    archived_best_A_.clear();
    CreateInitialPopulation();
    x_vectors_next_generation_ = x_vectors_current_;
    EvaluateCurrentVectors();
    evaluated_fitness_for_next_generation_ =
      evaluated_fitness_for_current_vectors_;
    for (long g = 0; g < total_generations_max_; ++g) {
      if (process_rank_ == kOutput && g%100 == 0) printf("%li\n",g);      
      to_be_archived_best_A_.clear();
      successful_mutation_parameters_S_F_.clear();
      successful_crossover_parameters_S_CR_.clear();        
      // //debug section
      // if (process_rank_ == kOutput)
      //   printf("==============  Generation %li =============\n", g);
      // PrintPopulation();      
      // PrintEvaluated();
      //end of debug section
      for (unsigned long i = 0; i < subpopulation_; ++i) {
        SetCRiFi(i);
        std::vector<double> mutated_v, crossover_u;
        mutated_v = Mutation(i);
        crossover_u = Crossover(mutated_v, i);
        Selection(crossover_u, i);
      }  // end of for all individuals in subpopulation
      ArchiveCleanUp();
      Adaption();
      x_vectors_current_.swap(x_vectors_next_generation_);
      evaluated_fitness_for_current_vectors_
        .swap(evaluated_fitness_for_next_generation_);
      SortEvaluatedCurrent();
      if (error_status_) return error_status_;
    }  // end of stepping generations
    PrintPopulation();      
    PrintEvaluated();
    return kDone;
  }  // end of int SubPopulation::RunOptimization()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> SubPopulation::Crossover(std::vector<double> mutation_v, long i) {
    const double CR_i = crossover_CR_[i];
    std::vector<double> crossover_u, x_current;
    crossover_u.resize(dimension_);
    x_current = x_vectors_current_.at(i);
    unsigned long j_rand = randint(0, dimension_ - 1);
    for (unsigned long c = 0; c < dimension_; ++c) {
      if (c == j_rand || rand(0,1) < CR_i)
        crossover_u[c] = mutation_v[c];
      else
        crossover_u[c] = x_current[c];
    }
    // //debug section
    // if (process_rank_ == kOutput)
    //   printf("x -> v -> u with CR_i=%4.2f j_rand=%li\n", CR_i, j_rand);
    // PrintSingleVector(mutation_v);
    // PrintSingleVector(x_current);
    // PrintSingleVector(crossover_u);
    // //end of debug section
    return crossover_u;
  } // end of  std::vector<double> SubPopulation::Crossover();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> SubPopulation::Mutation(long i) {
    std::vector<double> mutation_v, x_best_current, x_random_current,
      x_random_archive_and_current, x_current;    
    x_current = x_vectors_current_.at(i);
    x_best_current = GetXpBestCurrent();
    // //debug
    // if (process_rank_ == kOutput) printf("x_best: ");
    // PrintSingleVector(x_best_current);
    long index_of_random_current = -1;
    x_random_current = GetXRandomCurrent(&index_of_random_current, i);
    // //debug
    // if (process_rank_ == kOutput) printf("x_random: ");
    // PrintSingleVector(x_random_current);    
    x_random_archive_and_current = GetXRandomArchiveAndCurrent(index_of_random_current, i);
    // //debug
    // if (process_rank_ == kOutput) printf("x_random with archive: ");
    // PrintSingleVector(x_random_archive_and_current);
    mutation_v.resize(dimension_);
    double F_i = mutation_F_[i];
    for (unsigned long c = 0; c < dimension_; ++c) {
      // Mutation
      mutation_v[c] = x_current[c]
        + F_i * (x_best_current[c] - x_current[c])
        + F_i * (x_random_current[c]
                 - x_random_archive_and_current[c]);
      // Bounds control
      if (mutation_v[c] > x_ubound_[c])
        mutation_v[c] = (x_ubound_[c] + x_current[c])/2;
      if (mutation_v[c] < x_lbound_[c])
        mutation_v[c] = (x_lbound_[c] + x_current[c])/2;
    }
    // //debug section
    // int isSame = 888, isSame2 = 7777, isSame3 =11111;
    // for (long c = 0; c < dimension_; ++c) {
    //   double norm = std::abs(mutation_v[c] - x_current[c]);
    //   double norm2 = std::abs(x_random_current[c] - x_current[c]);
    //   double norm3 = std::abs(x_random_current[c]
    //                           - x_random_archive_and_current[c]);
    //   if ( norm  > 0.0001) isSame = 0;
    //   if ( norm2  > 0.0001) isSame2 = 0;
    //   if ( norm3  > 0.0001) isSame3 = 0;
    // }
    // if (process_rank_ == kOutput) printf("mutation_v%i%i%i:  ",isSame, isSame2, isSame3);
    // PrintSingleVector(mutation_v);
    // if (process_rank_ == kOutput) printf("current _v%i%i%i:  ",isSame, isSame2, isSame3);
    // PrintSingleVector(x_current);    
    // if (process_rank_ == kOutput)
    //   printf("  -> f = %4.2f                                    F_i=%4.2f\n",
    //          FitnessFunction(mutation_v), F_i);
    // //end of debug section
    return mutation_v;
  } // end of std::vector<double> SubPopulation::Mutation();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> SubPopulation::GetXpBestCurrent() {
    const unsigned long n_best_total = static_cast<long>
      (floor(subpopulation_ * best_share_p_ ));
    if (n_best_total == subpopulation_) error_status_ = kError; //TODO change kError to throw exception
    long best_n = randint(0, n_best_total);
    long best_n_index = -1, i = 0;
    for (auto x : evaluated_fitness_for_current_vectors_) {
      if (i == best_n) {
        best_n_index = x.second;
        break;
      }
      ++i;
    }
    return x_vectors_current_.at(best_n_index);
  }  // end of std::vector<double> SubPopulation::GetXpBestCurrent();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> SubPopulation::GetXRandomCurrent(long *index, long forbidden_index) {
    long random_n = randint(0, subpopulation_-1);
    while (random_n == forbidden_index) random_n = randint(0, subpopulation_-1);
    (*index) = random_n;
    return x_vectors_current_.at(random_n);
  }  // end of std::vector<double> SubPopulation::GetXRandomCurrent()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> SubPopulation::GetXRandomArchiveAndCurrent
  (unsigned long forbidden_index1, unsigned long forbidden_index2) {
    long archive_size = archived_best_A_.size();
    unsigned long random_n = randint(0, subpopulation_ + archive_size - 1);
    while (random_n == forbidden_index1 || random_n == forbidden_index2)
      random_n = randint(0, subpopulation_ + archive_size - 1);
    if (random_n < subpopulation_) return x_vectors_current_.at(random_n);
    random_n -= subpopulation_;
    unsigned long i = 0;
    for (auto x : archived_best_A_) {
      if (i == random_n) {
        // //debug
        // if (process_rank_ == kOutput) printf("Using Archive!!\n");
        return x;
      }
      ++i;
    }  // end of selecting from archive
    error_status_ = kError; //TODO change kError to throw exception
    std::vector<double> x;
    return x;    
  }  // end of std::vector<double> SubPopulation::GetXRandomArchiveAndCurrent()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::SetCRiFi(long i) {
    long k = 0;
    while (1) {
      mutation_F_[i] = randc(adaptor_mutation_mu_F_, 0.1);
      if (mutation_F_[i] > 1) {
        mutation_F_[i] = 1;
        break;
      }
      if (mutation_F_[i] > 0) break;
      ++k;
      if (k > 10) {
        mutation_F_[i] = 0.001;
        break;
      }
      if (k > 100) printf("k");
    }
    crossover_CR_[i] = randn(adaptor_crossover_mu_CR_,0.1);
    if (crossover_CR_[i] > 1) crossover_CR_[i] = 1;
    if (crossover_CR_[i] < 0) crossover_CR_[i] = 0;    
    return kDone;
  }  // end of int SubPopulation::SetCRiFi(long i)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// %todo Change returned kError to meangfull error code.
  int SubPopulation::Init(long total_population, long dimension) {// NOLINT
    total_population_  = total_population;
    if (total_population_ < 1)
      throw std::invalid_argument("You should set population > 0!");
    dimension_ = dimension;
    if (dimension_ < 1) 
      throw std::invalid_argument("You should set dimension > 0!");
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes_);
    if (process_rank_ < 0)
      throw std::invalid_argument("MPI problem: process_rank_ < 0!");
    if (number_of_processes_ < 1)
      throw std::invalid_argument("MPI problem: number_of_processes_ < 1!");
    std::random_device rd;
    generator_.seed(rd());

    // //debug
    // CheckRandom();

    subpopulation_ = total_population;

    current_generation_ = 0;
    x_vectors_current_.resize(subpopulation_);
    for (auto &x : x_vectors_current_) x.resize(dimension_);
    x_vectors_next_generation_.resize(subpopulation_);
    for (auto &x : x_vectors_next_generation_) x.resize(dimension_);
    evaluated_fitness_for_current_vectors_.resize(subpopulation_);
    // //debug
    // if (process_rank_ == kOutput) printf("%i, x1 size = %li \n", process_rank_, x_vectors_current_.size());
    x_lbound_.resize(dimension_);
    x_ubound_.resize(dimension_);
    mutation_F_.resize(subpopulation_);
    crossover_CR_.resize(subpopulation_);
    return kDone;
  }  // end of void SubPopulation::Init()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::PrintPopulation() {
    if (process_rank_ == kOutput) {
      printf("\n");	
      for (auto x : evaluated_fitness_for_current_vectors_) {
	long n = x.second;
	double fitness = x.first;
        printf("%6.2f:% 3li||", fitness, n);
        for (unsigned long c = 0; c < dimension_; ++c) 
          printf(" %+7.2f ", x_vectors_current_[n][c]);
	printf("\n");      
      }

      // for (long i = 0; i < subpopulation_; ++i) {
      //   printf("n%li:", i);
      //   for (long c = 0; c < dimension_; ++c) {
      //     printf(" %5.2f ", x_vectors_current_[i][c]);
      //   }
      //   printf("\n");
      // }  // end of for each individual
    }  // end of if output
    return kDone;
  }  // end of int SubPupulation::PrintPopulation()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::PrintEvaluated() {
    if (process_rank_ == kOutput) {
      for (auto x : evaluated_fitness_for_current_vectors_)
        printf("%li:%4.2f  ", x.second, x.first);
      printf("\n");
    }  // end of if output
    return kDone;
  }  // end of int SubPopulation::PrintEvaluated()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::PrintSingleVector(std::vector<double> x) {
    if (process_rank_ == kOutput) {
      for (auto c : x) printf("%5.2f ", c);
      printf("\n");
    }  // end of output
    return kDone;
  }  // end of int SubPopulation::PrintSingleVector(std::vector<double> x)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void SubPopulation::SetFeed(std::vector<std::vector<double> > x_feed_vectors) {
    isFeed_  = true;
    if (dimension_ < 1) 
      throw std::invalid_argument("You should set dimension before feed!");
    for (auto x : x_feed_vectors)
      if (x.size() != dimension_) 
	throw std::invalid_argument
	  ("Feed and optimization dimenstion should be the same size!");
    x_feed_vectors_.clear();
    for (auto &x : x_feed_vectors)
      x_feed_vectors_.push_back(x);
  }  // end of void SubPopulation::SetFeed()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::CreateInitialPopulation() {
    if ((subpopulation_ - x_feed_vectors_.size()) <1)
      throw std::invalid_argument("Too large feed!");
    x_vectors_current_.resize(subpopulation_ - x_feed_vectors_.size());
    for (auto &x : x_vectors_current_)
      for (unsigned long i = 0; i < dimension_; ++i) {
        if (x_lbound_[i] > x_ubound_[i])
	  throw std::invalid_argument("Wrong order of bounds!");
        x[i] = rand(x_lbound_[i], x_ubound_[i]);                            // NOLINT
      }  // end of for each dimension
    // //debug
    // for (auto x : x_vectors_current_[0]) if (process_rank_ == kOutput) printf("%g ",x);
    for (auto x: x_feed_vectors_) {
      x_vectors_current_.push_back(x);
        if (process_rank_ == kOutput) {
	  printf("--=-- Feed:\n");
	  for (auto index:x) printf(" %+7.2f", index);
	  printf("\n");
	}	
    }
    if (x_vectors_current_.size() != subpopulation_)
      throw std::invalid_argument("Population is not full after feed!");	
    x_feed_vectors_.clear();
    // if (process_rank_ == kOutput) {
    //   printf("==--== Initial population:\n");    
    //   EvaluateCurrentVectors();
    //   PrintPopulation();
    // }
    return kDone;
  }  // end of int SubPopulation::CreateInitialPopulation()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::SortEvaluatedCurrent() {
    evaluated_fitness_for_current_vectors_
      .sort([=](const std::pair<double, long>& a,                     // NOLINT
               const std::pair<double, long>& b) {                   // NOLINT
              bool cmp = a.first < b.first;
              if (is_find_minimum_) return cmp;
              return !cmp;
            });
    return kDone;
  }  // end of int SubPopulation::SortEvaluatedCurrent()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::EvaluateCurrentVectors() {
    evaluated_fitness_for_current_vectors_.clear();
    for (unsigned long i = 0; i < subpopulation_; ++i) {                        // NOLINT
      auto tmp = std::make_pair(FitnessFunction(x_vectors_current_[i]), i);
      evaluated_fitness_for_current_vectors_.push_back(tmp);
    }
    SortEvaluatedCurrent();
    // //debug
    // if (process_rank_ == kOutput) printf("\n After ");
    // for (auto val : evaluated_fitness_for_current_vectors_)
    //   if (process_rank_ == kOutput) printf("%g ", val.first);
    return kDone;
  }  // end of int SubPopulation::EvaluateCurrentVectors()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void SubPopulation::SetAllBoundsVectors
  (std::vector<double> lbound, std::vector<double> ubound) {
    x_lbound_.clear();
    x_ubound_.clear();    
    for (auto x : lbound) x_lbound_.push_back(x);
    for (auto x : ubound) x_ubound_.push_back(x);
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::SetAllBounds(double lbound, double ubound) {
    if (lbound >= ubound) {
      error_status_ = kError;
      return kError; //TODO change kError to throw exception
    }
    for (auto &x : x_lbound_) x = lbound;
    for (auto &x : x_ubound_) x = ubound;
    return kDone;
  }  // end of int SubPopulation::SetAllBounds(double lbound, double ubound)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double SubPopulation::randn(double mean, double stddev) {
    std::normal_distribution<double> distribution(mean, stddev);
    return distribution(generator_);
  }  // end of double SubPopulation::randn(double mean, double stddev)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double SubPopulation::randc(double location, double scale) {
    std::cauchy_distribution<double> distribution(location, scale);
    return distribution(generator_);
  }  // end of double SubPopulation::randc(double location, double scale)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  long SubPopulation::randint(long lbound, long ubound) {    // NOLINT
    std::uniform_int_distribution<long> distribution(lbound, ubound);  // NOLINT
    return distribution(generator_);
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double SubPopulation::rand(double lbound, double ubound) {                // NOLINT
    std::uniform_real_distribution<double> distribution(lbound, ubound);
    return distribution(generator_);
  }  // end of double rand(double lbound, double ubound)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void SubPopulation::CheckRandom() {
    std::map<int, int> hist_n, hist_c, hist_int, hist;
    int factor = 1000;
    for (int n = 0; n < 100*factor; ++n) {
      ++hist_n[std::round(randn(0, 2))];
      ++hist_c[std::round(randc(0, 2))];
      ++hist_int[std::round(randint(0, 10))];
      ++hist[std::round(rand(0, 15))];                                      // NOLINT 
    }
    if (process_rank_ == kOutput) {
      printf("Normal (0,2)\n");
      for (auto p : hist_n) {
        if (p.second > factor)
          printf("%i: % 4i %s\n", process_rank_, p.first,
                 std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist_n
      printf("Cauchy (0,2)\n");
      for (auto p : hist_c) {
        if (p.second > factor)
          printf("%i: % 4i %s\n", process_rank_, p.first,
                 std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist_c
      printf("Uniform int (0,10)\n");
      for (auto p : hist_int) {
        if (p.second > factor)
          printf("%i: % 4i %s\n", process_rank_, p.first,
                 std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist_int
      printf("Uniform double (0,15)\n");
      for (auto p : hist) {
        if (p.second > factor)
          printf("%i: % 4i %s\n", process_rank_, p.first,
                 std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist
    }  //  end of if current process_rank_ == kOutput
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::SetBestShareP(double p) {
    if (p < 0 || p > 1) {
      error_status_ = kError;
      return kError; //TODO change kError to throw exception
    }
    best_share_p_ = p;
    return kDone;
  }  // end of int SubPopulation::SetBestShareP(double p)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::SetAdapitonFrequencyC(double c) {
    if (c < 0 || c > 1) {
      error_status_ = kError;
      return kError; //TODO change kError to throw exception
    }
    adaptation_frequency_c_ = c;
    return kDone;
  }  // end of int SubPopulation::SetAdapitonFrequency(double c)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::SetDistributionLevel(int level) {
    distribution_level_ = level;
    return kDone;
  }
  // end of int SubPopulation::SetDistributionLevel(int level)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::PrintParameters(std::string comment) {
    if (process_rank_ == 0) {
      printf("#%s dim=%li NP=%li(of %li) p=%4.2f c=%4.2f generation=%li\n",
             comment.c_str(),dimension_, subpopulation_, total_population_,
             best_share_p_, adaptation_frequency_c_, total_generations_max_
             );
      fflush(stdout);
    }  // end of output
    return kDone;
  }  // end of int SubPopulation::PrintParameters(std::string comment)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::PrintResult(std::string comment) {    
    if (distribution_level_ == 0) {      
      auto x = evaluated_fitness_for_current_vectors_.front();
      std::vector<double> to_send {x.first};
      //printf("%8.5g\n  ", x.first);
      AllGatherVectorDouble(to_send);
      if (process_rank_ == 0) {
        double sum = 0;
        double size = static_cast<double>(recieve_double_.size());
        for (auto x : recieve_double_) sum += x;
        double mean = sum/size;
        double sigma = 0;
        for (auto x : recieve_double_) sigma += pow2(x - mean);
        sigma = sqrt(sigma/size);
        printf("%s gen%li, mean %4.1e (stddev %4.1e = %3.2g %%) runs(%g)\n",
               comment.c_str(), total_generations_max_,  mean,sigma,
               sigma*100.0/mean,size);
        // for (auto x : recieve_double_)
        //   printf("%18.15g\n", x);
      }
    }
    return kDone;
  }  // end of int SubPopulation::PrintResult()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> SubPopulation::GetFinalFitness() {
    recieve_double_.clear();    
    if (distribution_level_ == 0) {      
      auto x = evaluated_fitness_for_current_vectors_.front();
      std::vector<double> to_send {x.first};
      AllGatherVectorDouble(to_send);      
    }
    return recieve_double_;
  }  // end of sdt::vector<double> SubPopulation::GetFinalFitness();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> SubPopulation::GetBest(double *best_fitness) {
    auto x = evaluated_fitness_for_current_vectors_.front();
    (*best_fitness) = x.first;
    return x_vectors_current_[x.second];
  }  // end of std::vector<double> SubPopulation::GetBest(double *best_fitness)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> SubPopulation::GetWorst(double *worst_fitness) {
    auto x = evaluated_fitness_for_current_vectors_.back();
    (*worst_fitness) = x.first;
    return x_vectors_current_[x.second];
  }  // end of std::vector<double> SubPopulation::GetWorst(double *worst_fitness)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::AllGatherVectorDouble(std::vector<double> to_send) {
    long size_single = to_send.size();
    long size_all = size_single * number_of_processes_;
    recieve_double_.clear();
    recieve_double_.resize(size_all);
    MPI_Allgather(&to_send.front(), size_single, MPI_DOUBLE,
                  &recieve_double_.front(), size_single, MPI_DOUBLE,
                  MPI_COMM_WORLD);
    return kDone;
  }  // end of int SubPopulation::AllGatherVectorDouble(std::vector<double> to_send);
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::AllGatherVectorLong(std::vector<long> to_send) {
    long size_single = to_send.size();
    long size_all = size_single * number_of_processes_;
    recieve_long_.clear();
    recieve_long_.resize(size_all);
    MPI_Allgather(&to_send.front(), size_single, MPI_LONG,
                  &recieve_long_.front(), size_single, MPI_LONG,
                  MPI_COMM_WORLD);
    return kDone;
  }  // end of int SubPopulation::AllGatherVectorLong(std::vector<long> to_send);
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
}  // end of namespace jade
