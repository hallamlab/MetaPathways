#include "utilities.h"
//#include <iostream>
//#include <fstream>

std::string extract_sequence_name(const std::string &name) {
     char  cstr[1000000];
     strcpy(cstr, name.c_str());
     
     char * cptr = cstr;

     while( *cptr != '\t' && *cptr !=  ' ' && *cptr != '\0' )  cptr++; 
     (*cptr) ='\0';

     std::string sname(cstr);
     return sname;
}

void split(const string  &strn, std::vector<char *> &v, char *buf, char d) {
  strcpy(buf, strn.c_str());
  char *s1 = buf;
  v.clear();
  v.push_back(s1);
  while(*s1 != '\0') {
     if(*s1==d) { 
       *s1 = '\0';
       v.push_back(s1+1);
     }
     s1++;
  }
}

/*
std::string get_orf_name(std::string  &strn, std::vector<char *> &v, char *buf) {
    split(strn, v, buf, ';'); 
    if(v.size() == 0)  return std::string("");
    split(std::string(v[0]), v, buf, '=');
    if(v.size() < 2)  return std::string("");
    return std::string(v[1]);
}
*/

bool matchString(const string &str, const string & stringtomatch, bool fromstart) {
    unsigned long pos = str.find(stringtomatch);
    if(fromstart &&  pos ==0 ) return true;

    if( !fromstart && pos >= 0) return true;
    return false;
}

// Print system error and exit
void error(char *msg) {
    perror(msg);
    exit(1);
}

string shorten_id(string name, enum IDTYPE idtype)  {
  regex contigid_regex("_([0-9]+)$");
  regex orfid_regex("_([0-9]+_[0-9]+)$");
  regex orfid_trna_regex("_([0-9]+_[0-9]+_tRNA)$");
  regex orfid_rrna_regex("_([0-9]+_[0-9]+_rRNA)$");
  std::smatch sm; 
  switch (idtype) {
     case ORFID: 
       // prodigal predicted orf
       std::regex_search(name, sm, orfid_regex);
       if (sm.size() > 1) {
         return sm[1];
       }
       // tRNA orf
       std::regex_search(name, sm, orfid_trna_regex);
       if (sm.size() > 1) {
         return sm[1];
       }
       // rRNA orf
       std::regex_search(name, sm, orfid_rrna_regex);
       if (sm.size() > 1) {
         return sm[1];
       }
       break;
     case CONTIGID: 
       std::regex_search(name, sm, contigid_regex);
       break;
     default: 
       break;
  } 

  if (sm.size() > 1) {
     return sm[1];
  }
  return string("");
}

uint32_t strToUint16(string str1) {
  int myint(std::stoi(str1));
  uint32_t myint16(0);
  if (myint <= static_cast<int>(UINT16_MAX) && myint >=0) {
     myint16 = static_cast<uint32_t>(myint);
  } else {
     exception a;
     throw a;
  }
  return myint16;
}

void print_contig_orf_map(CONTIG_ORF *contig_orf_map) {
  for (auto it = contig_orf_map->begin(); it != contig_orf_map->end(); it++) {
    std::cout << it->first << std::endl;
    for (auto it2 = it->second->begin(); it2 != it->second->end(); it2++) {
      std::cout << "\t" <<  (*it2)->id << "\t" << (*it2)->start 
                << "\t" << (*it2)->end << "\t" << (*it2)->count << std::endl;
    }
  }
}

void sort_contig_orf_map(CONTIG_ORF *contig_orf_map) {
  for (auto it = contig_orf_map->begin(); it != contig_orf_map->end(); it++) {
    std::sort(it->second->begin(), it->second->end(), compare_orf_info); 
  }
}

void print_contig_read_counts(std::map<string, uint32_t> * contig_read_counts) {
  std::cout << "\nCONTIG\tREAD_COUNTS" << std::endl;
  for (auto it = contig_read_counts->begin(); it != contig_read_counts->end(); it++) {
    std::cout << it->first << "\t" << it->second << std::endl;
  }
}

bool compare_orf_info(ORFINFO *a, ORFINFO *b) { 
  return (a->start < b->start); 
}

bool compare_pairs_byid(const std::pair<string, float> &a, 
                        const std::pair<string, float> &b) { 
  return (a.first < b.first); 
}

void print_estimates(std::ostream *output, 
                     vector<std::pair<string, float>> &counts) {
  for (auto it1 = counts.begin(); it1 != counts.end(); it1++) {  
       *output << std::fixed << std::setprecision(2) 
               <<  it1->first << "\t" << it1->second <<  std::endl;
  }
}

void print_all_estimates(std::ostream *output, 
                     vector<std::pair<string, float>> &counts1,
                     vector<std::pair<string, float>> &counts2,
                     vector<std::pair<string, float>> &counts3) {
  for (auto it1 = counts1.begin(), it2 = counts2.begin(), it3= counts3.begin(); 
       it1 != counts1.end() && it2 != counts2.end() && it3 != counts3.end(); 
       it1++, it2++, it3++) {  

       *output << std::fixed << std::setprecision(2) 
               <<  it1->first << "\t" << it1->second 
               <<  "\t" << it2->second << "\t" 
               << it3->second <<  std::endl;
  }
}
