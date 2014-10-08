#ifndef COMPUTE_COMMUTATORS_H_
#define COMPUTE_COMMUTATORS_H_

#include "compute-commutators-util.h"

#include <set>
#include <string>
#include <vector>

namespace compute_commutators {

typedef std::vector<int> term;
typedef compute_commutators_util::single_coeffs single_coeffs;
typedef compute_commutators_util::TermsToCoeffsMap TermsToCoeffsMap;

class ComputeCommutators {
 public:
  ComputeCommutators(int n, bool verbose);
  // Add initial terms in Hamiltonian (pq and pqrs) along with symbolic
  // representation of coefficients to the terms to coefficients maps with 
  // terms normal ordered (so swapping takes place as necessary), and then
  // remove complex conjugates from map.
  void AddInitialTerms();
  // A helper for AddInitialTerms that returns the coefficient corresponding to
  // the initial term, accounting for symmetries that will be tracked in
  // unique_coeffs set.
  // The coefficient should have the same form as the initial term (i.e., same
  // list of indices [p,-q] or [p,q,-r,-s]) since we have not done any swaps
  // or multiplications yet at this point.
  std::vector<single_coeffs> GetInitialSumCoeffs(
      std::vector<int> curr_coeff_term);
  // Helper for GetInitalSumCoeffs that checks if coefficient term has a
  // symmetrical permutation already seen before. Returns corresponding pointer.
  std::set<term>::iterator InitialCoeffSeen(const int& one, const int& two,
      const int& three, const int& four);
  // Prepare the order of terms for the final calculation of Trotter error.
  void InterleaveTerms();
  // Helper for InterleaveTerms to add terms to interleaved_order
  void AddTermToInterleavedOrder(const term& curr_term);
  // Calculate the Trotter error using the terms arranged in the order from
  // InterleaveTerms, and with the complex conjugates added back in.
  void CalculateTrotterError();
  // Helper for CalculateTrotterError to return term, its conjugate, its coeff.
  std::pair<std::vector<term>, std::vector<single_coeffs> > GetTermForTrotter(
      const int& index);
  void PrintFinalResults(FILE* output);
 private:
  int num_orbitals;
  TermsToCoeffsMap initial_terms_to_coefficients;
  std::set<term> initial_terms;
  std::set<term> unique_coeffs;
  std::vector<term> interleaved_order;
  TermsToCoeffsMap final_terms_to_coefficients;
  bool verbose;
};

} //  namespace compute_commutators

#endif
