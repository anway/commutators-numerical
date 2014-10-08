#ifndef COMPUTE_COMMUTATORS_UTIL_H_
#define COMPUTE_COMMUTATORS_UTIL_H_

#include <map>
#include <set>
#include <string>
#include <vector>

namespace compute_commutators_util {

// Each term is represented by a list of indices, where a positive index
// represents a raising operator and a negative index represents a lowering
// operator.
typedef std::vector<int> term;

// single_coeff represents a term of the form integer_number * h_1 * h_2 * ...
// where the h_i are symbolic representations of h_{pq}, h_{pqrs}
// We will then represent a sum of single_coeff with a vector of single_coeffs
struct single_coeffs {
  int integer_multiplier = 1;
  std::multiset<std::vector<int> > product_of_coeffs;
};

class ComputeCommutatorsUtil {
 public:
  // Returns the conjugate of a term in normal order, or an empty list is a term
  // is its own conjugate.
  static term GetConjugate(const term& curr_term);
  // Print a vector of ints
  static void PrintIndices(FILE* output, const term& curr_term);
  // Print a vector of single_coeffs
  static void PrintSumOfCoeffs(FILE* output,
      const std::vector<single_coeffs>& sum_of_coeffs);
  // Multiply two coefficient vectors term by term to get a new coefficient
  // vector.
  static std::vector<single_coeffs> MultiplySumOfCoeffs(
      const std::vector<single_coeffs>& first,
      const std::vector<single_coeffs>& second); 
  // Check if double commutator [A, [B, C]] is trivially zero
  static bool TriviallyCommutes(const term& first_term, const term& second_term,
      const term& third_term);
  // Helper to concatenate three terms (e.g., ABC)
  static term ConcatenateThreeTerms(const term& first_term,
      const term& second_term, const term& third_term);
};

// Wrapper class for terms_to_coeffs map
class TermsToCoeffsMap {
 public:
  // Adds a term and coefficient to the map, with the term in normal order
  // after necessary swapping.
  void AddNormalForm(term curr_term, std::vector<single_coeffs> curr_coeff);
  // Removes complex conjugates from the map.
  void RemoveComplexConjugates();
  bool HasTerm(const term& curr_term);
  // Return iterator to start of map.
  std::map<term, std::vector<single_coeffs> >::iterator Begin();
  // Return iterator to end of map.
  std::map<term, std::vector<single_coeffs> >::iterator End();
  // Returns value corresponding to key, or throws an exception if no value.
  std::vector<single_coeffs> At(const term& key);
 private:
  std::map<term, std::vector<single_coeffs> > terms_to_coefficients;
};

}  // namespace compute_commutators_util

#endif
