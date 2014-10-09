#include "compute-commutators.h"

int main(int argc, char* argv[])
{
  if (argc < 3 || argc > 4) {
    printf("Usage: ./compute-commutators molecule basis (verbose version) or\n"
        "./compute-commutators molecule basis (non-verbose version).\n");
    return(1);
  } 
  bool verbose;
  if (argc == 3) {
    verbose = true;
  } else if (*argv[3] == 'F' || *argv[3] == 'f') {
    verbose = false;
  } else {
    printf("Usage: ./compute-commutators molecule basis (verbose version) or\n"
        "./compute-commutators molecule basis (non-verbose version).\n");
    return(2);
  }

  std::string in_file_name = "data/from_jarrod/" + std::string(argv[1]) + "-" +
      std::string(argv[2]) + ".int";
  compute_commutators::ComputeCommutators compute_commutators =
      compute_commutators::ComputeCommutators(verbose);
  compute_commutators.AddInitialTerms(in_file_name);
  compute_commutators.InterleaveTerms();
  compute_commutators.CalculateTrotterError();

  // Print results.
  char out_file_name[80]; 
  sprintf(out_file_name, "data/error_terms/%s_%s.txt", argv[1], argv[2]);
  FILE *file_p = fopen(out_file_name, "w");
  compute_commutators.PrintFinalResults(file_p);
  return 0;
}
