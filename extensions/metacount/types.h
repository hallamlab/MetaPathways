#ifndef __DATA_TYPES__
#define __DATA_TYPES__

#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <functional>
#include <numeric>

using namespace std;

typedef struct _ORFINFO {
  std::string id;
  uint32_t start, end;
  float count;
  
  _ORFINFO(string _id, uint32_t _start, 
           uint32_t _end): id(_id), start(_start), 
           end(_end), count(0) {}
} ORFINFO;

enum IDTYPE {ORFID, CONTIGID};

typedef map<std::string, vector<ORFINFO *> *>  CONTIG_ORF;

typedef struct _TRIPLET {
    unsigned int start, end;
    float multi;
} TRIPLET;

typedef vector<TRIPLET> TRIPLETS;

typedef struct _CONTIG {
    unsigned int L;
    TRIPLETS M;
} CONTIG;

enum STATS_TYPE {TPM, RPKM, COUNT};
enum READNO {FIRST, SECOND};
enum ALNVALID {VALID, INVALID};
enum READTYPE {SE, PE};
enum ALIGNTYPE {NO_ALIGN, PRIMARY, NON_PRIMARY, OTHER};
enum CHIMERIC{SUPPLEMENTARY, NON_SUPPLEMENTARY};

typedef  struct _MATCH {
    std::string query, subject;
    int16_t start, end, alnlength;

    int16_t samflag;

    enum ALNVALID alnvalid;   // VALID, INVALID
    enum READNO readno;   //  FIRST, SECOND
    enum READTYPE readtype;   //;  SE, PE (single-end or paired-end)    
    enum ALIGNTYPE aligntype;  // PRIMARY, NON_PRIMARY
    enum CHIMERIC chimeric;  // SUPPLEMENTARY, NON_SUPPLEMENTARY
 
    float weight;

    _MATCH(): weight(0) { } 

    string print();
} MATCH;


typedef struct GLOBAL_STATS {
    int num_files; 
    int num_contigs; 
    int total_contig_length; 
    int total_reads1; 
    int total_reads2; 
    float tot_count_across_orfs;
    float tot_reads_across_orfs;

    GLOBAL_STATS();

    void print(std::ostream *output);
} GLOBAL_STATS;

typedef struct SAM_STATS {
    string file_name;
    int num_unmapped_alns ; 
    int num_mapped_alns ; 
    int num_reads_1 ; 
    int num_reads_2 ; 
    int num_reads1_on_orfs ; 
    int num_reads2_on_orfs ; 
    int num_total_alns ;
    
    SAM_STATS();

    struct SAM_STATS&  operator+( const struct SAM_STATS & stats);

    void print(std::ostream *output);
} SAM_STATS;

typedef struct _RESULTS {
   vector<SAM_STATS *> sam_file_results;
   GLOBAL_STATS global_stats;
   CONTIG_ORF *contig_orf;
   std::map<string, uint32_t> *contig_lengths;
   std::map<string, uint32_t> * contig_read_counts; 

   void compute_stats();
   void print(std::ostream *output);

} RESULTS;

#endif //__DATA_TYPES__



