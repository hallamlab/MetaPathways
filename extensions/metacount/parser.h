#ifndef __PARSER_H__
#define __PARSER_H__

#include <string>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "utilities.h"
#include "types.h"

using namespace std;

class MatchOutputParser {
  public:
    std::string filename;
    std::string format;
    std::ifstream input;
    char buf[1000];
    vector<char *> tempv;
    vector<char *> fields;

    virtual ~MatchOutputParser() = 0;

    MatchOutputParser(const std::string &filename, const std::string &format);

    unsigned long get_Num_Unmapped_Reads();

    virtual bool nextline(MATCH &match)=0;
  protected:
    unsigned long num_unmapped_reads;

};

//subclass of the MathOutputParser
class GffFileParser: virtual public MatchOutputParser {
  public:
    GffFileParser(const std::string &filename, const std::string &format);
    virtual bool nextline(MATCH &match);
    ~GffFileParser();
};

class ParserFactory {
  public:
    static MatchOutputParser * createParser(const std::string &filename, const std::string &format) {
      if(format.find("GFF") != string::npos ) {
         return new GffFileParser( filename, format);
      }
      return 0;
    }
};


void getSamFlagInfo(uint32_t c, MATCH &match);
#endif // __PARSER_H__
