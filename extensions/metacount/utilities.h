#ifndef __METACOUNT_UTILS__H__
#define __METACOUNT_UTILS__H__

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <cstdlib>
#include <regex>
#include <exception>

#include "types.h"

using namespace std;


void split(const std::string  &strn, std::vector<char *> &v, char *buf, char d='\t');

//std::string get_orf_name(std::string & strn, std::vector<char *> &v, char *buf);


bool matchString(const string &str, const string & stringtomatch, bool fromstart=false);

void get_fasta_sequence_info(const std::string &fasta_file_name);

std::string extract_sequence_name(const std::string &name);

string to_string(unsigned long i);

// Print system error and exit
void error(char *msg);

string shorten_id(string name, enum IDTYPE);

uint32_t strToUint16(string str1);

void print_contig_orf_map(CONTIG_ORF *);

void sort_contig_orf_map(CONTIG_ORF *contig_orf_map);

void print_contig_read_counts(std::map<string, uint32_t> * contig_read_counts);

void print_estimates(std::ostream *output, vector<std::pair<string, float>> &counts);

void print_all_estimates(std::ostream *output, 
                     vector<std::pair<string, float>> &counts1,
                     vector<std::pair<string, float>> &counts2,
                     vector<std::pair<string, float>> &counts3);

bool compare_pairs_byid(const std::pair<string, float> &a, 
                        const std::pair<string, float> &b);

bool compare_orf_info(ORFINFO *a, ORFINFO *b);  

#endif // __METACOUNT_UTILS__H__

