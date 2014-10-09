#include "compute-commutators-util.h"

#include <regex>
#include <set>

namespace compute_commutators_util {

term ComputeCommutatorsUtil::GetConjugate(const term& curr_term) {
  term conjugate;
  for (auto rit = curr_term.rbegin(); rit != curr_term.rend(); ++rit) {
    conjugate.push_back(-1 * (*rit));
  }
  if (conjugate == curr_term) {
    // Return nothing if a term is equal to its own conjugate.
    conjugate.clear();
  } 
  return conjugate;
}

void ComputeCommutatorsUtil::PrintIndices(FILE* output, const term& curr_term) {
  fprintf(output, "[ ");
  for (int index : curr_term) {
    fprintf(output, "%d ", index);
  }
  fprintf(output, "] ");
}

bool ComputeCommutatorsUtil::TriviallyCommutes(const term& first_term,
    const term& second_term, const term& third_term) {
  // First check if B and C have indices in common
  std::set<int> set_one(second_term.begin(), second_term.end());
  std::set<int> set_two(third_term.begin(), third_term.end());
  std::vector<int> intersection;
  std::set_intersection(set_one.begin(), set_one.end(), set_two.begin(),
      set_two.end(), std::back_inserter(intersection));
  if (!intersection.empty()) {
    return false;
  }
  intersection.clear();
  // Now check if A and BC have indices in common
  std::copy(third_term.begin(), third_term.end(), std::inserter(
      set_one, set_one.end()));
  std::set<int> set_three(first_term.begin(), first_term.end());
  std::set_intersection(set_one.begin(), set_one.end(), set_three.begin(),
      set_three.end(), std::back_inserter(intersection));
  if (!intersection.empty()) {
    return false;
  }
  return true;
}

term ComputeCommutatorsUtil::ConcatenateThreeTerms(const term& first_term,
    const term& second_term, const term& third_term) {
  term result;
  result.insert(result.end(), first_term.begin(), first_term.end());
  result.insert(result.end(), second_term.begin(), second_term.end());
  result.insert(result.end(), third_term.begin(), third_term.end());
  return result;
}

std::vector<std::string> ComputeCommutatorsUtil::Split(const std::string& str) {
  std::regex rgx("\s|\t");
  std::sregex_token_iterator
      first{begin(str), end(str), rgx, -1},
      last;

  return{first, last};
}

void TermsToCoeffsMap::AddNormalForm(term curr_term, double curr_coeff) {
  for (std::vector<int>::iterator it = curr_term.begin(); it != curr_term.end();
      ++it) {
    for (std::vector<int>::iterator rit = it; rit != curr_term.begin();
        --rit) {
      const int right = *rit, left = *(rit - 1);
      if (right > left) {
        // Swap the two indices.
        *(rit - 1) = right;
        *rit = left;
        if (left == -1 * right) {
          // If additionally we are swapping two indices that act on the same
          // tensor factor, add an additional term without the two indices.
          // (fermionic commutation relations)
          term exchange_term;
          exchange_term.insert(exchange_term.end(), curr_term.begin(), rit - 1);
          exchange_term.insert(exchange_term.end(), rit + 1, curr_term.end());
          if (!exchange_term.empty()) {
            AddNormalForm(exchange_term, curr_coeff);
          }
        }
        // Flip the sign of the coefficient because we swapped.
        curr_coeff *= -1;
      } else {  // No more terms out of order, so stop swapping.
        break;
      }
    }
  }
  std::set<int> no_repeats(curr_term.begin(), curr_term.end());
  // If an index is repeated, we have two raising or two lowering operators, so
  // don't add to map.
  if (no_repeats.size() == curr_term.size()) {
    if (terms_to_coefficients.find(curr_term) == terms_to_coefficients.end()) {
      // Term not in map already.
      terms_to_coefficients[curr_term] = curr_coeff;
    } else {  // Term already in map; just add (sum) existing coefficients.
      terms_to_coefficients[curr_term] += curr_coeff;
    }
  }
}

void TermsToCoeffsMap::RemoveComplexConjugates() {
  for (const auto& term_to_coeff : terms_to_coefficients) {
    term conjugate = compute_commutators_util::ComputeCommutatorsUtil
        ::GetConjugate(term_to_coeff.first);
    if (!conjugate.empty()) {  // If the term is not its own conjugate
      auto it = terms_to_coefficients.find(conjugate);
      // Remove its conjugate from map.
      if (it != terms_to_coefficients.end()) {
        terms_to_coefficients.erase(it);
      }
    }    
  }
}

bool TermsToCoeffsMap::HasTerm(const term& curr_term) {
  if (terms_to_coefficients.find(curr_term) != terms_to_coefficients.end()) {
    return true;
  } else {
    return false;
  }
}

std::map<term, double>::iterator TermsToCoeffsMap
    ::Begin() {
  return terms_to_coefficients.begin();
}

std::map<term, double>::iterator TermsToCoeffsMap
    ::End() {
  return terms_to_coefficients.end();
}

double TermsToCoeffsMap::At(const term& key) {
  return terms_to_coefficients.at(key);
}

}  // namespace compute_commutators_util
