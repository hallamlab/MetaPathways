#include "parser.h"

using namespace std;

MatchOutputParser::MatchOutputParser(const std::string &filename, const std::string &format) {
  this->filename = filename;
  this->format = format;
  this->num_unmapped_reads =0;
};

unsigned long MatchOutputParser::get_Num_Unmapped_Reads() {
  return  this->num_unmapped_reads;
}

MatchOutputParser::~MatchOutputParser() {
}

/**
    @brief getSamFlagInfo gets the samflag into
    @detail
    get thes the samflag info based on the following bit positions on a
    16 bit unsigned int (uint32_t)

    description of bit positions
        1 = read paired
        2 = read mapped in proper pair
        3 = read unmapped
        4 = mate unmapped
        5 = read reverse strand
        6 = mate reverse strand
        7 = first in pair
        8 = second in pair
        9 = not primary alignment
        10 = read fails platform/vendor quality checks
        11 = read is PCR or optical duplicate
        12 = supplementary alignment 

    Weights of the mapping are computed according to the following rule

    Mapping/alignment information:
    readtype: SE, PE (single-end or paired-end)    
    weight: 0, 0.5, 1
    aligntype: PRIMARY, NON_PRIMARY

    Notations:
    1st UNMAP: means read itself is not mapped
    2nd UNMAP: means the read's mate is not mapped

    Weights according to combinations <value>(<bit position>)
    a) SE + NO_ALIGN              0(1) + 4(3) = 4               : 0
    b) SE + PRIMARY               0(1) + 0(3) + 0(9)  = 0       : 1    

    c) PE + PROPER + UNMAP + MAP + PRIMARY  1(1) + 2(2) + 4(3) + 0(4) + 0(9) = 7         : 1
    d) PE + PROPER + MAP + UNMAP + PRIMARY  1(1) + 2(2) + 0(3) + 8(4) + 0(9) = 11        : 1

    e) PE + UNMAP + UNMAP   1(1) + 4(3) + 8(4)             : 0
    f) PE + PROPER + MAP + MAP + PRIMARY  1(1) + 2(2) + 0(3) + 0(4) + 0(9)     : 0.5

    g) rest 

    @params     c  flag 
    @params match  returns the information about the alignment
*/

void getSamFlagInfo(uint32_t c, MATCH &match)  {
    match.samflag = c;

    // supplementary
    if (((c >> 11) & 1) == 0) {
      match.chimeric = NON_SUPPLEMENTARY;
    } else {
      match.chimeric = SUPPLEMENTARY;
    }

    // se is FIRST and for pe it is dependent on the flag
    if ((c & 1) == 0) {
      match.readno = FIRST;
    } else {
      if ((c >> 6 & 1) == 1) {
         match.readno = FIRST;
      } else {
         match.readno = SECOND;
      }
    }

    match.alnvalid = INVALID;
    if (match.alnlength > 0) {
       match.alnvalid = VALID;
    }
    //a)
    if ( 
        ((c & 1) == 0) && 
        (((c >> 2) & 1) == 1) 
        ) {  // 0
        match.readtype = SE;
        match.aligntype = NO_ALIGN;
        match.weight = 0.0;
    // b)
    } else if (
        ((c & 1) == 0) && 
        (((c >> 2) & 1) == 0) 
        ) {  // 0
       match.readtype = SE;
       match.aligntype = PRIMARY;
       match.weight = 1.0;
    // c)
    } else if (
        ((c & 1) == 1) &&
        (((c >> 1) & 1) == 0) && 
        (((c >> 2) & 1) == 1) && 
        (((c >> 3) & 1) == 0) 
        ) {  // 0
      match.readtype = PE;
      if (match.alnlength > 0) {
         match.aligntype = PRIMARY;
         match.weight = 1.0;  //1.0
      } else {
         match.aligntype = NO_ALIGN;
         match.weight = 0.0;  //1.0
      }
    // d)
    } else if (
       ((c & 1) == 1) && 
       (((c >> 1) & 1) == 0) && 
       (((c >> 2) & 1) == 0) && 
       (((c >> 3) & 1) == 1) 
       ) {  // 0
       match.readtype = PE;
       if (match.alnlength > 0) {
          match.aligntype = PRIMARY;
          match.weight = 1.0;  //1.0
       } else {
           match.aligntype = NO_ALIGN;
           match.weight = 0.0;  //1.0
       }
    // e)
    } else if (
       ((c & 1) == 1) && 
       (((c >> 2) & 1) == 1) && 
       (((c >> 3) & 1) == 1) 
       ) {  // 0
       match.readtype = PE;
       match.aligntype = NO_ALIGN;
       match.weight = 0.0;   // 0;
    // f)
    } else if ((c & 1) == 1 && 
       ((c >> 1) & 1) == 1 && 
       ((c >> 2) & 1) == 0 && 
       ((c >> 3) & 1) == 0 
       ) {  // 0
        match.readtype = PE;
        match.aligntype = PRIMARY;
        match.weight = 0.5;  // 0.5;
    } else { // g)
        if ( (c & 1) == 1) {
            match.readtype = PE;
         } else { 
             match.readtype = SE;
         }
         match.aligntype = OTHER;
         match.weight = 0.0;
    }
}

GffFileParser::GffFileParser(const std::string &filename, 
                             const std::string &format) : MatchOutputParser(filename, format) {

  this->input.open(filename.c_str(), std::ifstream::in);
  if (!this->input.good()) {
    std::cerr << "Error opening '"<<filename<<"'. Bailing out." << std::endl;
    return ;
  }  
}

GffFileParser::~GffFileParser() {
   this->input.close();
}

/*

 regex orfname_long_regex("^ID=[^;]+_([0-9]+_[0-9]+);");
 regex orfname_short_regex("^ID=([0-9]+_[0-9]+);");
 regex orfname_short_trna_regex("^ID=([0-9]+_[0-9]+_tRNA);");
 regex orfname_short_rrna_regex("^ID=([0-9]+_[0-9]+_rRNA);");
 std::smatch sm; 
*/

namespace parser {
  std::smatch sm; 
 regex orfname_long_regex("^ID=[^;]+_([0-9]+_[0-9]+);");
 regex orfname_short_regex("^ID=([0-9]+_[0-9]+);");
 regex orfname_short_trna_regex("^ID=([0-9]+_[0-9]+_tRNA);");
 regex orfname_short_rrna_regex("^ID=([0-9]+_[0-9]+_rRNA);");
}
string get_orf_name(string name)  {
 // regex orfname_long_regex("^ID=[^;]+_([0-9]+_[0-9]+);");
  //regex orfname_short_regex("^ID=([0-9]+_[0-9]+);");
  //regex orfname_short_trna_regex("^ID=([0-9]+_[0-9]+_tRNA);");
 // regex orfname_short_rrna_regex("^ID=([0-9]+_[0-9]+_rRNA);");
//  std::smatch sm; 

  // a long name prodigal orf
  std::regex_search(name, parser::sm, parser::orfname_long_regex);
  if (parser::sm.size() > 1) {
     return parser::sm[1];
  }

  // a short name prodigal orf
  std::regex_search(name, parser::sm, parser::orfname_short_regex);
  if (parser::sm.size() > 1) {
     return parser::sm[1];
  }

  // a short name tRNA orf
  std::regex_search(name, parser::sm, parser::orfname_short_trna_regex);
  if (parser::sm.size() > 1) {
     return parser::sm[1];
  }

  // a short name rRNA orf
  std::regex_search(name, parser::sm, parser::orfname_short_rrna_regex);
  if (parser::sm.size() > 1) {
     return parser::sm[1];
  }

  return string("");
}


bool GffFileParser::nextline(MATCH &match) {
  string line;
  string field8; 
  bool _success = false;
  static int i = 0;
  
  while (std::getline(this->input, line ).good()) {
      split(line, fields, this->buf, '\t');
      //std::cout << line << std::endl;
     // std::cout << "size " << fields.size() << std::endl;
      if(fields.size() < 9)  continue;
      //std::cout << fields[3] << " - " << fields[4] << " " << _success <<std::endl;
      _success = true;
      break;
  }  
  i++;
  
  if (_success)  {  
     match.start = atoi(fields[3]);
     match.end = atoi(fields[4]);
     match.query = fields[0]; 

/*        if(false &&  i%100 ==0  ) {
        std::cout << "i= " <<  i << std::endl;
        std::cout << "subject " <<  match.subject << std::endl;
     }
*/
     //respect the order before calling get_orf_name
     // because the same buffer is used
     //field8 = std::string(fields[8]);
     //match.subject = get_orf_name(field8, this->tempv, this->buf);
     match.subject =  get_orf_name(std::string(fields[8]));
     return true;
  }
  
  return false;
}
