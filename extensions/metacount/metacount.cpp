#include <map>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <experimental/filesystem>
#include <functional>
#include <numeric>


#include "types.h"
#include "utilities.h"
#include "options.h"
#include "core.h"
#include "SamValidation.h"

using namespace std;
namespace fs = std::experimental::filesystem;


bool compare_triplets(const TRIPLET &a, const TRIPLET &b) {
  return a.start < b.start ? true : false; 
}
 
int main(int argc, char **argv ){
  // parse options
  Options options;

  std::ostream *err = &std::cerr;
  // if (argc < 9) { options.print_usage(argv[0]); exit(1); }
  if (options.SetOptions(argc, argv)==false) { 
    options.print_usage(argv[0]); exit(0); 
  }

  RESULTS results;
  

  *err << "Creating contig to orf map....";
  results.contig_orf = create_contig_orf_map(options.orf_file);

  for (size_t i = 0; i < options.read_map_files.size(); i++) {
     if (!fs::exists(options.read_map_files[i].c_str())) {
        std::cout << "ERROR: file" << options.read_map_files[i] << " is missing!\n";
        std::cerr << "ERROR: File " << options.read_map_files[i] << " is missing!\n";
        exit(0);
     }
  }
  *err << "done" << std::endl;

  results.global_stats.num_contigs = results.contig_orf->size();
  results.global_stats.num_files = options.read_map_files.size();

  *err << "Reading contig information....";
  results.contig_lengths = get_contig_information(options.read_map_files);
  *err << "done" << std::endl;

  //std::cout << "num contigs " << contig_lengths->size() << std::endl;

  //print_contig_orf_map(results.contig_orf);

  *err << "Sorting contig to orf map....";
  sort_contig_orf_map(results.contig_orf);
  *err << "done" << std::endl;

  *err << "Reading bam files...";
  results.contig_read_counts = 
    read_bam_files(results, options.read_map_files);
  *err << "done" << std::endl;

  // compute the stats
  *err << "Computing stats...";
  results.compute_stats();
  *err << "done" << std::endl;


  std::ostream *output = &std::cout;
  std::ofstream fout;
  if (options.output_file.size() > 0) {
      try {
         fout.open(options.output_file.c_str(), std::ifstream::out);
      }
      catch(exception& e) {
         error(e.what());
      }
      output = &fout;
  }
  
  vector<std::pair<string, float>> orf_count, orf_count_rpkm, orf_count_tpm;
  if (options.count_type == string("TPM")) {
     orf_count = compute_stats(results.contig_orf, results, TPM);
     *output << "#ORFID\tTPM" << std::endl;
     print_estimates(output, orf_count);
  } else if (options.count_type == string("COUNT")) {
     orf_count = compute_stats(results.contig_orf, results, COUNT);
     *output << "#ORFID\tCOUNT" << std::endl;
     print_estimates(output, orf_count);
  } else if (options.count_type == string("RPKM")) {
     orf_count = compute_stats(results.contig_orf, results, RPKM);
     *output << "#ORFID\tRPKM" << std::endl;
     print_estimates(output, orf_count);
  } else if (options.count_type == string("ALL")) {
     orf_count = compute_stats(results.contig_orf, results, COUNT);
     orf_count_tpm = compute_stats(results.contig_orf, results, RPKM);
     orf_count_rpkm = compute_stats(results.contig_orf, results, TPM);
     *output << "#ORFID\tCOUNT\tTPM\tRPKM" << std::endl;
     print_all_estimates(output, orf_count, orf_count_tpm, orf_count_rpkm);
  } else {
     error("ERROR: unknown type of estimates suggested");
  }

  if (options.show_stats) { 
      std::ostream *output = &std::cout;
      std::ofstream fout1;
      if (options.stats_file.size() > 0) {
         try {
            fout1.open(options.stats_file.c_str(), std::ifstream::out);
         }
         catch(exception& e) {
            error(e.what());
         }
         output = &fout1;
      }   
      results.print(output);
      fout1.close();
  }

  /*
  print_contig_read_counts(results.contig_read_counts);
  print_contig_orf_map(results.contig_orf);
*/

  exit(0);
   




  return 0;
} 
