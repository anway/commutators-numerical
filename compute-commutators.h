#ifndef COMPUTE_COMMUTATORS_H_
#define COMPUTE_COMMUTATORS_H_

#include "compute-commutators-util.h"

#include <set>
#include <string>
#include <vector>

namespace compute_commutators {

typedef std::vector<int> term;
typedef compute_commutators_util::TermsToCoeffsMap TermsToCoeffsMap;

class ComputeCommutators {
 public:
  ComputeCommutators(bool verbose);
  // Add initial terms in Hamiltonian (pq and pqrs) along with numeric
  // coefficients to the terms to coefficients maps with 
  // terms normal ordered (so swapping takes place as necessary), and then
  // remove complex conjugates from map. Takes terms and coefficients from input
  // file formatted as index index coeff, or index index index index coeff
  void AddInitialTerms(std::string file_name);
  // Prepare the order of terms for the final calculation of Trotter error.
  void InterleaveTerms();
  // Helper for InterleaveTerms to add terms to interleaved_order
  void AddTermToInterleavedOrder(const term& curr_term);
  // Calculate the Trotter error using the terms arranged in the order from
  // InterleaveTerms, and with the complex conjugates added back in.
  void CalculateTrotterError();
  // Helper for CalculateTrotterError to return term, its conjugate, its coeff.
  std::pair<std::vector<term>, double> GetTermForTrotter(const int& index);
  void PrintFinalResults(FILE* output);
 private:
  int num_orbitals = 0;
  TermsToCoeffsMap initial_terms_to_coefficients;
  std::set<term> initial_terms;
  std::set<term> unique_coeffs;
  std::vector<term> interleaved_order;
  TermsToCoeffsMap final_terms_to_coefficients;
  bool verbose;
};

} //  namespace compute_commutators

#endif
