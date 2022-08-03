#include "types.h"

string MATCH::print() {
  string out(this->query);
  out = out + "\t";
  out = out + std::to_string(this->samflag);
  out = out + "\t";
  out = out + this->subject;
  out = out + "\t";
  out = out + std::to_string(this->start);
  out = out + "\t";
  out = out + std::to_string(this->end);
  out = out + "\t";
  out = out + std::to_string(this->alnlength);
  out = out + "\t";
  if (this->readtype == SE) {
      out = out + "SE";
  } else { 
      out = out + "PE";
  }
  out = out + "\t";
  if (this->readno == FIRST) {
      out = out + "R1";
  } else { 
      out = out + "R2";
  }
  out = out + "\t";
  if (this->alnvalid == VALID) {
      out = out + "VALID";
  } else { 
      out = out + "INVALID";
  }
  out = out + "\t";
  switch (this->aligntype) {
      case NO_ALIGN:
          out = out + "\t" + "NO_ALIGN";
          break;
      case PRIMARY:
          out = out + "\t" + "PRIMARY";
          break;
      case NON_PRIMARY:
          out = out + "\t" + "NON_PRIMARY";
          break;
      default:
          out = out + "\t" + "OTHER";
          break;
  } 

  switch (this->chimeric) {
      case SUPPLEMENTARY:
          out = out + "\t" + "SUPPLEMENTARY";
          break;
      case NON_SUPPLEMENTARY:
          out = out + "\t" + "NON_SUPPLEMENTARY";
          break;
      default:
          out = out + "\t" + "OTHER";
  }

  out = out + "\t" +  std::to_string(this->weight);
  return out;
} 


SAM_STATS::SAM_STATS() {
        num_unmapped_alns =0; 
        num_mapped_alns =0; 
        num_reads_1 =0; 
        num_reads_2 =0; 
        num_total_alns = 0;
        num_reads1_on_orfs = 0 ; 
        num_reads2_on_orfs = 0; 
    }

struct SAM_STATS& SAM_STATS::operator+( const struct SAM_STATS & stats) {
        this->num_unmapped_alns += stats.num_unmapped_alns;  
        this->num_mapped_alns  += stats.num_mapped_alns  ;
        num_reads_1 += stats.num_reads_1;  
        num_reads_2 += stats.num_reads_2; 
        num_total_alns  += stats.num_total_alns  ;
        return *this;
    }


void SAM_STATS:: print(std::ostream *output) {
        *output << std::endl;
        *output << "Name of file                     : " << file_name << std::endl; 
        *output << "Number of alignments             : " << num_total_alns << std::endl; 

        int x = num_total_alns == 0 ? 1 : num_total_alns;           

        *output << "Number of unmapped alignments    : " << num_unmapped_alns << std::endl;
        *output << "Number of mapped alignments      : " << num_mapped_alns << std::endl; 

        *output << "Percentage of mapped alignments  : " << std::setprecision(4)  
                << static_cast<float>(num_mapped_alns) *100/static_cast<float>(x) << "%" <<  std::endl; 
        *output << "Number of read 1                 : " << num_reads_1 << std::endl; 
        *output << "Number of read 2                 : " << num_reads_2 << std::endl; 
        *output << "Number of read1 on ORFs          : " << num_reads1_on_orfs << std::endl; 
        *output << "Number of read2 on ORFs          : " << num_reads2_on_orfs << std::endl; 
}

GLOBAL_STATS::GLOBAL_STATS() {
      num_files = 0; 
      num_contigs = 0; 
      total_contig_length = 0; 
      total_reads1 = 0; 
      total_reads2 = 0; 
      tot_count_across_orfs = 0;
}


void RESULTS::compute_stats() {
  global_stats.total_contig_length = 
       accumulate(
                   contig_lengths->begin(), 
                   contig_lengths->end(), 0, 

                   [] (uint32_t value, const std::pair<string, uint32_t>&p) {
                        return value + p.second; 
                   }
                 );

  global_stats.tot_count_across_orfs
           = accumulate(
                 contig_orf->begin(), 
                 contig_orf->end(), 0, 
                    
                 [] (uint32_t value, const std::pair<string, vector<ORFINFO *>*> &orfvec) {
                       float tot_count = value;
                       for (auto it = orfvec.second->begin(); it != orfvec.second->end(); it++) {
                          tot_count = tot_count + (*it)->count; 
                       }
                       return tot_count; 
                 }
               );

  for (auto it = this->sam_file_results.begin(); it != sam_file_results.end(); it++) {
      global_stats.total_reads1 += (*it)->num_reads_1;
      global_stats.total_reads2 += (*it)->num_reads_2;
  }

}

void GLOBAL_STATS::print(std::ostream *output) {
        *output << std::endl;
        *output << "Stats across all SAM files" <<  std::endl; 
        *output << "Number of sam/bam files           : " << num_files << std::endl; 
        *output << "Number of contigs                 : " << num_contigs << std::endl; 
        *output << "Total contig length               : " << total_contig_length << std::endl; 
        *output << "Total read1s                      : " << total_reads1 << std::endl;
        *output << "Total read2s                      : " << total_reads2 << std::endl;
        *output << "Total count to across orfs        : " << tot_count_across_orfs << std::endl; 
}

void RESULTS::print(std::ostream *output) {
    global_stats.print(output);

    for (auto it = this->sam_file_results.begin(); it != sam_file_results.end(); it++) {
      (*it)->print(output);
    }

}
