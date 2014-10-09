#ifndef COMPUTE_COMMUTATORS_UTIL_H_
#define COMPUTE_COMMUTATORS_UTIL_H_

#include <map>
#include <string>
#include <vector>

namespace compute_commutators_util {

// Each term is represented by a list of indices, where a positive index
// represents a raising operator and a negative index represents a lowering
// operator.
typedef std::vector<int> term;

class ComputeCommutatorsUtil {
 public:
  // Returns the conjugate of a term in normal order, or an empty list is a term
  // is its own conjugate.
  static term GetConjugate(const term& curr_term);
  // Print a vector of ints
  static void PrintIndices(FILE* output, const term& curr_term);
  // Check if double commutator [A, [B, C]] is trivially zero
  static bool TriviallyCommutes(const term& first_term, const term& second_term,
      const term& third_term);
  // Helper to concatenate three terms (e.g., ABC)
  static term ConcatenateThreeTerms(const term& first_term,
      const term& second_term, const term& third_term);
  // Helper to split a string on a space delimiter.
  static std::vector<std::string> Split(const std::string& str);
};

// Wrapper class for terms_to_coeffs map
class TermsToCoeffsMap {
 public:
  // Adds a term and coefficient to the map, with the term in normal order
  // after necessary swapping.
  void AddNormalForm(term curr_term, double curr_coeff);
  // Removes complex conjugates from the map.
  void RemoveComplexConjugates();
  bool HasTerm(const term& curr_term);
  // Return iterator to start of map.
  std::map<term, double>::iterator Begin();
  // Return iterator to end of map.
  std::map<term, double>::iterator End();
  // Returns value corresponding to key, or throws an exception if no value.
  double At(const term& key);
 private:
  std::map<term, double> terms_to_coefficients;
};

}  // namespace compute_commutators_util

#endif
