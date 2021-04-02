#include "options.h"
//#include <iostream>
//#include <fstream>


void Options::print_usage(char *arg) {
  std::cout << "NAME\n"   
             << "\t" << "metacount - computes the TPM/RPKM stats for ORF from  aligned sam/bam files" \
             "\n\n"
             <<"SYNOPSIS\n " 
             << "\t" << "metacount --sam [SAM/BAM] [--sam [SAM/BAM]] -gff [GFF] [--print-stats] " 
             " [--stats-out-file]\n" "\t[--stats-out-file FILE] [--estimate-type [TYPE]] [--out-file FILE]\n\n" \
             <<"DESCRIPTION\n "  
             << "\t" << "computes the TPM/RPKM stats for ORF from  aligned sam/bam files\n" \
             << "\t-h\n\t\tprint help\n"\
             << "\t--sam\n\t\te.g. --samsam/bam [--sam sam/bam]  [at least one REQUIRED] \n"\
             << "\t--gff\n\t\tGFF file with orf name in the last column as \"ID=samplename_42_0\" [REQUIRED]\n"
             << "\t--print-stats\n\t\tprints overall stats [OPTIONAL]\n"\
             << "\t--stats-out-file\n\t\tfile name to where stats are expected if --print-stats [OPTIONAL]\n"\
             << "\t--estimate-type\n\t\tTYPE: TPM/RPKM/COUNT/ALL [OPTIONAL default TPM]\n"\
             << "\t--out-file\n\t\toutput file name  [OPTIONAL default stdout] \n"\
             << std::endl;
}

bool Options::SetOptions(int argc, char *argv[]) { 
    for(int i = 1; i < argc ; i++) {   
        if (strncmp(argv[i], "-h", strlen("-h")) == 0) {
            print_usage(argv[0]);
            exit(0);
        }   
        else if (strncmp(argv[i], "--estimate-type", strlen("--estimate-type")) == 0) {
            this->count_type = argv[++i];
            if (this->count_type != "ALL" && this->count_type != "TPM" 
                && this->count_type != "RPKM" && this->count_type != "COUNT") {
               print_usage(argv[0]);
               exit(0);
            }
        }   
        else if (strncmp(argv[i], "--stats-out-file", strlen("--stats-out-file")) == 0) {
           if (i+1 < argc) {
             this->stats_file = argv[++i];
           } else {
             print_usage(argv[0]);
             exit(0);
           }
        }   
        else if (strncmp(argv[i], "--sam", strlen("--sam")) == 0) {
            this->read_map_files.push_back(string(argv[++i]));
        }   
        else if (strncmp(argv[i], "--print-stats", strlen("--print-stats")) == 0) {
            this->show_stats = true;
        }   
        else if ( strncmp(argv[i], "--out-file", strlen("--out-file")) == 0) {
            this->output_file = argv[++i];
        } 
        else if (strncmp(argv[i], "--gff", strlen("--gff")) == 0) {
            this->orf_file = argv[++i];
        }
        else {
            std::cout << "ERROR: Input options : " << argv[i] <<  std::endl;
            print_usage(argv[0]);
            exit(0);
        }
    } //for loop for arguments processing

    bool exit_with_error = false;
    if (this->read_map_files.size() == 0) {
        std::cout << "ERROR: There must be a at least SAM/BAM file" << std::endl;
        exit_with_error = true;
    }
    
    if (this->orf_file.size()==0) {
        std::cout << "ERROR: There must be an ORF info file [GFF format]" << std::endl;
        exit_with_error = true;
    }
    if (exit_with_error) {
       std::cout <<  "Type \"metacount -h\" for more info\n";
       exit(1);
    }

    return true;
};

void Options::print_options() {
    std::cout << "Alignment file (SAM format)  : "<< this->read_map_files.size() << std::endl; 
    std::cout << "ORF file [GFF]                    : "<< this->orf_file << std::endl; 
    std::cout << "Output file                       : "<< this->output_file << std::endl; 
}

