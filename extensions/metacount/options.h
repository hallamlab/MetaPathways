#ifndef __OPTIONS__H__
#define __OPTIONS__H__
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <cstdlib>
#include <regex>

#include "types.h"

using namespace std;

// Structure for RPKM input options
struct Options {
    /* Input files for RPKM */
    string orf_file; // .gff ORF file from PRODIGAL
    string output_file; // location to write output file (i.e., update pathway table)
    string stats_file; // location to write output file (i.e., update pathway table)
    bool show_stats;
    vector<string> read_map_files;
    
    /* Flags and settings */
    bool read_counts; // show counts not rpkm
    string reads_map_file_format; // aligner type BWA or BLAST, two SAM files or one
    
    string count_type;
    Options(){
       stats_file = ""; // location to write output file (i.e., update pathway table)
       read_map_files.clear();
       orf_file = "";
       output_file = "";

       show_stats = false;
       read_counts = false;
       string count_type = "TPM";
    };
    
    void print_usage( char *arg);
    void print_options();
 
    bool SetOptions(int argc, char *argv[]);
	// bool SetOptionCommon( const char *flag, const char *value );
	// bool SetOption( const char *flag, const char *value );
	// bool SetOption2D( const char *flag, const char *value );
	// bool SetOptionEST( const char *flag, const char *value );
	// bool SetOptions( int argc, char *argv[], bool twodata=false, bool est=false );

	void Validate();
	// void ComputeTableLimits( int naan, int typical_len, size_t mem_need );

	void Print();
};
#endif // __OPTIONS__H__

