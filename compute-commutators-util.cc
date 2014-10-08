#include "compute-commutators-util.h"

namespace compute_commutators_util {

typedef compute_commutators_util::single_coeffs single_coeffs;

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

void ComputeCommutatorsUtil::PrintSumOfCoeffs(FILE* output,
    const std::vector<single_coeffs>& sum_of_coeffs) {
  for (auto it = sum_of_coeffs.begin(); it != sum_of_coeffs.end();
      ++it) {
    fprintf(output, "%d * ", (*it).integer_multiplier);
    for (const auto& prod_of_terms : (*it).product_of_coeffs) {
      PrintIndices(output, prod_of_terms);
    }  
    if (it != sum_of_coeffs.end() - 1) {
      fprintf(output, " + ");
    }
  }
}

std::vector<single_coeffs> ComputeCommutatorsUtil::MultiplySumOfCoeffs(
    const std::vector<single_coeffs>& first,
    const std::vector<single_coeffs>& second) {
  // result will carry result of multiplication
  std::vector<single_coeffs> result;
  for (const single_coeffs& first_term : first) {
    for (const single_coeffs& second_term : second) {
      // first multiply just the symbolic terms (without coefficient)
      std::multiset<std::vector<int> > new_term;
      new_term.insert(first_term.product_of_coeffs.begin(),
          first_term.product_of_coeffs.end());
      new_term.insert(second_term.product_of_coeffs.begin(),
          second_term.product_of_coeffs.end());    
      // check if resulting symbolic term has already been seen in this
      // multiplication
      std::vector<single_coeffs>::iterator it_result;
      for (it_result = result.begin(); it_result != result.end();
          ++it_result) {
        if (it_result->product_of_coeffs == new_term) {
          break;
        }
      }
      // if not, push new symbolic term and coefficient
      if (it_result == result.end()) {
        single_coeffs new_term_result;
        new_term_result.integer_multiplier = first_term.integer_multiplier *
            second_term.integer_multiplier;
        new_term_result.product_of_coeffs = new_term;
        result.push_back(new_term_result);
      } else {  // if so, just add to coefficient of existing term
        it_result->integer_multiplier += first_term.integer_multiplier *
            second_term.integer_multiplier; 
      }
    }
  }
  return result;
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
  } else {
    return true;
  }
}

term ComputeCommutatorsUtil::ConcatenateThreeTerms(const term& first_term,
    const term& second_term, const term& third_term) {
  term result;
  result.insert(result.end(), first_term.begin(), first_term.end());
  result.insert(result.end(), second_term.begin(), second_term.end());
  result.insert(result.end(), third_term.begin(), third_term.end());
  return result;
}

void TermsToCoeffsMap::AddNormalForm(term curr_term,
    std::vector<single_coeffs> curr_coeff) {
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
        for (auto it = curr_coeff.begin(); it != curr_coeff.end(); ++it) {
          it->integer_multiplier *= -1;
        }
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
      for (const auto& current_coeff : curr_coeff) {
        if (current_coeff.integer_multiplier != 0) {
          std::vector<single_coeffs>::iterator it;
          // Check if we're just combining existing terms or adding new term.
          for (it = terms_to_coefficients[curr_term].begin();
              it != terms_to_coefficients[curr_term].end(); ++it) {
            if (it->product_of_coeffs == current_coeff.product_of_coeffs) {
              break;
            }
          }
          // If not, push new term and coefficient.
          if (it == terms_to_coefficients[curr_term].end()) {
            terms_to_coefficients[curr_term].push_back(current_coeff);
          } else {
            it->integer_multiplier += current_coeff.integer_multiplier;
            if (it->integer_multiplier == 0) {
              terms_to_coefficients[curr_term].erase(it);
            }
          }
          }
      }
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

std::map<term, std::vector<single_coeffs> >::iterator TermsToCoeffsMap
    ::Begin() {
  return terms_to_coefficients.begin();
}

std::map<term, std::vector<single_coeffs> >::iterator TermsToCoeffsMap
    ::End() {
  return terms_to_coefficients.end();
}

std::vector<single_coeffs> TermsToCoeffsMap::At(const term& key) {
  return terms_to_coefficients.at(key);
}

}  // namespace compute_commutators_util
