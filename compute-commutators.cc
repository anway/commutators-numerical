#include "compute-commutators.h"
#include "compute-commutators-util.h"

#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <regex>

namespace compute_commutators {

ComputeCommutators::ComputeCommutators(bool verbose) : verbose(verbose) {} 

void ComputeCommutators::AddInitialTerms(std::string file_name) {
  std::string line;
  std::ifstream file(file_name);
  // Check for lines that are formatted as indices and coefficient.
  std::regex regex_filter("\[0-9]+\.\[0-9]+|\[0-9]+");
  if (file.is_open()) {
    while (std::getline(file, line)) {
      if (std::regex_search(line, regex_filter)) {
        std::vector<std::string> split_line = compute_commutators_util
            ::ComputeCommutatorsUtil::Split(line);
        if (split_line.size() == 3) {
          term curr_term;
          // Note that orbitals are 0-indexed so we need to add 1 to
          // everything
          const int& p = std::stoi(split_line[0], nullptr, 10) + 1;
          const int& q = std::stoi(split_line[1], nullptr, 10) + 1;
          const double& curr_coeff = std::stof(split_line[2], nullptr);
          curr_term.push_back(p);
          curr_term.push_back(-1 * q);

          initial_terms_to_coefficients.AddNormalForm(curr_term, curr_coeff);

          num_orbitals = std::max(num_orbitals, p);
          num_orbitals = std::max(num_orbitals, q);
        } else if (split_line.size() == 5) {
          term curr_term;
          // Note that orbitals are 0-indexed so we need to add 1 to
          // everything
          const int& p = std::stoi(split_line[0], nullptr, 10) + 1;
          const int& q = std::stoi(split_line[1], nullptr, 10) + 1;
          const int& r = std::stoi(split_line[2], nullptr, 10) + 1;
          const int& s = std::stoi(split_line[3], nullptr, 10) + 1;
          const double& curr_coeff = std::stof(split_line[4], nullptr) / -2.0;
          curr_term.push_back(p);
          curr_term.push_back(q);
          curr_term.push_back(-1 * r);
          curr_term.push_back(-1 * s);

          initial_terms_to_coefficients.AddNormalForm(curr_term, curr_coeff);

          num_orbitals = std::max(num_orbitals, p);
          num_orbitals = std::max(num_orbitals, q);
          num_orbitals = std::max(num_orbitals, r);
          num_orbitals = std::max(num_orbitals, s);
        } 
      }
    }
    file.close();
  } else {
    fprintf(stderr, "Unable to open file.\n");
  }

  initial_terms_to_coefficients.RemoveComplexConjugates();

  // Put all the keys into initial_terms.
  for (auto it = initial_terms_to_coefficients.Begin();
      it != initial_terms_to_coefficients.End(); ++it) {
    initial_terms.insert(it->first);
  }
}

void ComputeCommutators::AddTermToInterleavedOrder(const term& curr_term) {
  if (initial_terms_to_coefficients.HasTerm(curr_term)) {
    auto it = initial_terms.find(curr_term);
    if (it != initial_terms.end()) { 
      interleaved_order.push_back(curr_term);
      initial_terms.erase(initial_terms.find(curr_term));
    }
  }      
}

void ComputeCommutators::InterleaveTerms() {
  // First add Hpp terms.
  for (int p = 1; p <= num_orbitals; ++p) {
    term curr_term;
    curr_term.push_back(p);
    curr_term.push_back(-1 * p);

    AddTermToInterleavedOrder(curr_term);
  } 
  // Now add Hpqqp terms.
  for (int p = 1; p <= num_orbitals; ++p) {
    for (int q = 1; q <= p; ++q) {
      term curr_term;
      curr_term.push_back(p);
      curr_term.push_back(q);
      curr_term.push_back(-1 * q);
      curr_term.push_back(-1 * p);

      AddTermToInterleavedOrder(curr_term);
    }
  }
  // Now interleave (p, -q) and (p, r, -r, -q) terms.
  for (int p = 1; p <= num_orbitals; ++p) {
    for (int q = 1; q <= num_orbitals; ++q) {
      term curr_term;
      curr_term.push_back(p);
      curr_term.push_back(-1 * q);
  
      AddTermToInterleavedOrder(curr_term);

      for (int r = 1; r <= p && r <= q; ++r) {
        term curr_term;
        curr_term.push_back(p);
        curr_term.push_back(r);
        curr_term.push_back(-1 * r);
        curr_term.push_back(-1 * q); 

        AddTermToInterleavedOrder(curr_term);
      }
    }
  }
  // Now add pqrs terms.
  for (int p = 1; p <= num_orbitals; ++p) {
    for (int q = 1; q <= p; ++q) {
      for (int r = 1; r <= num_orbitals; ++r) {
        for (int s = r; s <= num_orbitals; ++s) {
          term curr_term;
          curr_term.push_back(p);
          curr_term.push_back(q);
          curr_term.push_back(-1 * r);
          curr_term.push_back(-1 * s);
           
          AddTermToInterleavedOrder(curr_term);
        }
      }
    }
  }

  if (verbose) {
    fprintf(stderr, "Order of terms in Trotter series:\n");
    for (const auto& curr_term : interleaved_order) {
      fprintf(stderr, "Term: ");
      compute_commutators_util::ComputeCommutatorsUtil::PrintIndices(stderr,
          curr_term);
      fprintf(stderr, " Coeff: %f\n",
          initial_terms_to_coefficients.At(curr_term));
    }
  }
  fprintf(stderr, "Number of terms in interleaved order: %lu\n",
      interleaved_order.size());
}

std::pair<std::vector<term>, double> ComputeCommutators
    ::GetTermForTrotter(const int& index) {
  term curr_term = interleaved_order[index];
  term curr_conjugate = compute_commutators_util::ComputeCommutatorsUtil
      ::GetConjugate(curr_term);
  double curr_coeff = initial_terms_to_coefficients.At(curr_term);
  std::vector<term> term_and_conjugate;
  if (!curr_conjugate.empty()) {
    // return both term and conjugate if conjugate and term differ
    term_and_conjugate.push_back(curr_term);
    term_and_conjugate.push_back(curr_conjugate);
  } else {
    // otherwise return just term
    term_and_conjugate.push_back(curr_term);
  }
  return std::make_pair(term_and_conjugate, curr_coeff);
}

void ComputeCommutators::CalculateTrotterError() {
  int num_commutators = 0, num_terms = interleaved_order.size();
  for (int i = 0; i < num_terms; ++i) {
    num_commutators += i * (i + 1);  
  }
  int one_percent = num_commutators / 100.0;
  fprintf(stderr, "There are %d possible commutators.\n", num_commutators);
  int counter = 0;
  clock_t start = clock();

  // Loop over the possible combinations a <= b, c < b of
  // (1 / 12) * [A * (1 - delta(A, B)/2), [B, C]] = ...
  // ((1 - delta(A, B) / 2) / 12) * (ABC - ACB - BCA + CBA)

  // Loop over B.
  for (int b = 1; b < num_terms; ++b) {
    const auto& B_term_coeff = GetTermForTrotter(b);
    const std::vector<term>& B_terms = B_term_coeff.first;
    const double&  B_coeff = B_term_coeff.second;

    // Loop over A.
    for (int a = 0; a <= b; ++a) {
      const auto& A_term_coeff = GetTermForTrotter(a);
      const std::vector<term>& A_terms = A_term_coeff.first;
      const double& A_coeff = A_term_coeff.second;

      // Loop over C.
      for (int c = 0; c < b; ++c) {
        const auto& C_term_coeff = GetTermForTrotter(c);
        const std::vector<term>& C_terms = C_term_coeff.first;
        const double& C_coeff = C_term_coeff.second;
       
        // Increment counter.
        counter += 1;

        // Calculate the coefficient (multiplying the 3 terms). 
        double multiplied_coeffs = A_coeff * B_coeff * C_coeff;
        // To handle the scale (each coefficient needs to be multipled by 1 / 24
        // if a == b, or 1 / 12 if not): we multiply the coefficient by 2 if
        // a is not equal to b, and divide everything by 24 in the very end.
        // This allows us to keep coefficients as ints. 
        if (a == b) {
          multiplied_coeffs /= 24.0; 
        } else {
          multiplied_coeffs /= 12.0;
        }

        // Compute commutators.
        for (term A : A_terms) {
          for (term B : B_terms) {
            for (term C : C_terms) {
              if (!compute_commutators_util::ComputeCommutatorsUtil
                  ::TriviallyCommutes(A, B, C)) {
                // ABC term
                final_terms_to_coefficients.AddNormalForm(
                    compute_commutators_util::ComputeCommutatorsUtil::
                    ConcatenateThreeTerms(A, B, C), multiplied_coeffs);
                // ACB term
                final_terms_to_coefficients.AddNormalForm(
                    compute_commutators_util::ComputeCommutatorsUtil::
                    ConcatenateThreeTerms(A, C, B), -1 * multiplied_coeffs);
                // BCA term
                final_terms_to_coefficients.AddNormalForm(
                    compute_commutators_util::ComputeCommutatorsUtil::
                    ConcatenateThreeTerms(B, C, A), -1 * multiplied_coeffs);
                // CBA term
                final_terms_to_coefficients.AddNormalForm(
                    compute_commutators_util::ComputeCommutatorsUtil::
                    ConcatenateThreeTerms(C, B, A), multiplied_coeffs);
              }
            }
          }
        }
      
        // Report progress.
        if (one_percent != 0 && counter % one_percent == 0) {
          int percent_complete = counter / one_percent;
          int elapsed = double(clock() - start) / CLOCKS_PER_SEC; 
          int rate = elapsed / percent_complete; 
          int eta = rate * (100 - percent_complete); 
          // Get current date and time.
          time_t rawtime;
          time(&rawtime);
          struct tm *timeinfo = localtime(&rawtime);
          fprintf(stderr, "%sComputation %d%% complete. Approximately "
              "%d minute(s) remaining.\n", asctime(timeinfo),
              percent_complete, eta / 60);
        }
      } 
    }
  } 
}

void ComputeCommutators::PrintFinalResults(FILE* output) {
  for (auto it = final_terms_to_coefficients.Begin();
      it != final_terms_to_coefficients.End(); ++it) {
    compute_commutators_util::ComputeCommutatorsUtil::PrintIndices(output,
        it->first);
    fprintf(output, " %.16e\n", it->second);
  }
}

}  // namespace compute-commutators
