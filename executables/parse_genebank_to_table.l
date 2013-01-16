%{
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define MAX_BUF 400000
#define MIN_BUF 100000
#define SMALL_BUF 1000
#define MAX_GENES 40000
#define PRINT 0


/* define all your global variable here */
int locus_count = 0;
char LOCUS_str[MAX_BUF];
char DEFINITION_str[MAX_BUF];
char ACCESSION_str[MAX_BUF];
char VERSION_str[MAX_BUF];
char KEYWORDS_str[MAX_BUF];
char SOURCE_str[MAX_BUF];
char ORGANISM_str[MAX_BUF];
char REFERENCE_str[MAX_BUF];
char AUTHORS_str[MAX_BUF];
char CONSRTM_str[MAX_BUF];
char TITLE_str[MAX_BUF];
char JOURNAL_str[MAX_BUF];
char PUBMED_str[MAX_BUF];
char COMMENT_str[MAX_BUF];
char FEATURES_str[MAX_BUF];
char L_source_str[MAX_BUF];
char F_strain_str[MAX_BUF];
char F_organism_str[MAX_BUF];
char F_chromosome_str[MAX_BUF];
char gene_str[MAX_BUF];
char F_codon_start_str[MAX_BUF];
char F_locus_tag_str[MAX_BUF];
char CDS_str[MAX_BUF];
char F_product_str[MAX_BUF];
char F_codon_str[MAX_BUF];
char F_transl_table_str[MAX_BUF];
char F_note_str[MAX_BUF];
char F_db_xref_str[MAX_BUF];
char F_EC_number_str[MAX_BUF];
char F_translation_str[MAX_BUF];
char ORIGIN_str[MAX_BUF];
char devNull_str[MAX_BUF];


/*define the functions */
void print_values();
void print_value( const char *, char *);
void print_value_delim(char *, char *);
void reset_values();
void partial_reset_values();
void clean_values();
void remove_digits(char *);
void remove_characters(char *,const char *);
void collapse_multiple_spaces_to_single(char *string);
void format_into_plain_string(char *string) ;
unsigned int first_locus_flag = 1;

/* this is a structure that hold the information about one gene/CDS */
struct __GENE {
  int a;
  char gene_str[MIN_BUF];
  char codon_start[MIN_BUF];
  char locus_tag[MIN_BUF];
  char CDS_str[MIN_BUF];
  char product[MIN_BUF];
  char transl_table[MIN_BUF];
  char note[MIN_BUF];
  char db_xref[MIN_BUF];
  char EC_number[MIN_BUF];
  char translation[MIN_BUF];
  struct __GENE *next;
};

typedef struct __GENE * GENEPTR;
typedef struct __GENE  GENE;

/* a maxinum of 20 genes are expected for one LOCUS/ORIGIN, typically it will
 rarely be above 2 */
GENE genes[MAX_GENES];

/* keeps track of the number of genes in one locus */
int gene_count = -1;


char * get_field_number(char *string, int i ) {
   char *word;
   char temp_string[MAX_BUF];;
   int j=1;
   strcpy(temp_string, string);
   word=strtok(temp_string, " "); 
   while( word !=NULL) {
     if(i==j) return(word);
     j++;
     word=strtok(NULL, " "); 
   }
}

void print_labels() {
     print_value_delim("LOCUS ID", "\t");
     print_value_delim("LOCUS SIZE","\t");
     print_value_delim("locus_tag","\t");
     print_value_delim("gene location","\t"); 
     print_value_delim("product","\t");
     print_value_delim("EC_number","\t");
     print_value_delim("note","\t");
     print_value_delim("db_xref","\n");
}

/* print alll the values realated to one locus */
void print_values(){
   char *word;

   /*print the genes/CDS if there are (0 or many ) */
  int i;

  /* first clean up the strings */
 // printf("Gene count = %d\n",gene_count);
  for(i=0; i<=gene_count; i++) {
     format_into_plain_string(genes[i].gene_str);
     format_into_plain_string(genes[i].CDS_str);
     format_into_plain_string(genes[i].product);
     format_into_plain_string(genes[i].codon_start);
     format_into_plain_string(genes[i].locus_tag);
     format_into_plain_string(genes[i].transl_table);
     format_into_plain_string(genes[i].EC_number);
     format_into_plain_string(genes[i].note);
     format_into_plain_string(genes[i].db_xref);

    /* use remove character function to remove even single spaces, new lines,
    tabs and digits */ 
     remove_characters(genes[i].translation,"\n\t 0123456789/\"");
  }
  /* print the values in the form of a table */
/*
      print_value_delim("LOCUS ID", "\t");
      print_value_delim("LOCUS SIZE","\t");
      print_value_delim("gene location","\t"); 
      print_value_delim("product","\t");
      print_value_delim("EC_number","\t");
      print_value_delim("note","\t");
      print_value_delim("db_xref","\n");
*/

  for(i=0; i<=gene_count; i++) {
      print_value_delim(get_field_number(LOCUS_str,1), "\t");
      print_value_delim(get_field_number(LOCUS_str,2),"\t");
      print_value_delim(genes[i].locus_tag,"\t");
      print_value_delim(genes[i].gene_str,"\t"); 
      print_value_delim(genes[i].product,"\t");
      print_value_delim(genes[i].EC_number,"\t");
      print_value_delim(genes[i].note,"\t");
      print_value_delim(genes[i].db_xref,"\n");
  }
}


/* format into a string with no new lines and multiple long spaces, 
   if there are multiple spaces then collapse them into a single space */
void format_into_plain_string(char *string) {
      remove_characters(string,"\n\t\"");
      collapse_multiple_spaces_to_single(string);
}

void clean_values() {
   remove_characters(F_translation_str,"\n\t 0123456789/\"");
   //remove_digits(ORIGIN_str);
   remove_characters(ORIGIN_str,"\n\t 0123456789/\"");
}

/* remove the characters in "characters" from the  the string called "string" for all such 
  occurances */
void remove_characters(char *string, const  char *characters) {
    int i,j,k;
    unsigned int flag;
    i=0; j=0; k=0;
    while(string[i]!='\0') {
      k=0;
      flag =0;
      while( characters[k]!='\0' ) {
         if(string[i] == characters[k] ) {
          flag=1;break;
         }
         k++;
      }
      if(flag==0) { 
        string[j]=string[i];j++;
      }
      i++;
    }
    string[j]='\0';
}

/*   if there are multiple spaces then collapse them into a single space */
void collapse_multiple_spaces_to_single(char *string){
    int i,j;
    unsigned space_pos=0;
    i=0;
    j=0;
    while(string[i]!='\0') {
      if(string[i] == ' ') 
          space_pos++;
      else
          space_pos=0;

      if(space_pos < 2  ) {
          string[j]= string[i];
          j++; 
      }
      i++;
    }
    string[j]='\0';
}

/* removes the digits from a string */
void remove_digits( char *string) {
    int i,j;
    i=0;
    j=0;
    while(string[i]!='\0') {
      if(string[i] < '0' ||  string[i] >'9' ) {
          string[j]= string[i];
          j++; 
      }
      i++;
    }
    string[j]='\0';
}

char *str_ptr;
void print_and_reset() {
    clean_values();
    print_values(); 
    reset_values(); 
}

%}

/* regular definitions */
delim    [ \t\n]
ws       {delim}+
letter   [A-Za-z]
digit    [0-9]
character    [A-Za-z0-9]  
sentence    [A-Za-z0-9\,\.]  
begspace ^[ ]*
space  [ \t]
transgram   [A-Za-z0-9:\- \n\"]
anything [.\n ]
anybutquotes [^"]

%s LOCUS
%s DEFINITION
%s ACCESSION
%s VERSION
%s KEYWORDS
%s SOURCE
%s ORGANISM
%s REFERENCE
%s AUTHORS
%s CONSRTM
%s TITLE
%s JOURNAL
%s PUBMED
%s COMMENT
%s FEATURES
%s L_source
%s F_strain
%s F_organism
%s F_chromosome
%s gene
%s F_codon_start
%s gene_F_locus_tag
%s F_locus_tag
%s CDS
%s F_product
%s F_codon
%s F_transl_table
%s F_note
%s F_db_xref
%s F_EC_number
%s F_translation
%s ORIGIN
%s  NULL_STATE
%%



"LOCUS"         {            
                        if( PRINT )  printf("%s\t",yytext);
                        if(first_locus_flag==0) {
                            clean_values();
                 /*           fprintf(yyout,"\n\n===================================================================\n"); */
                            print_values(); 
                            reset_values(); 
                        }
                        if(first_locus_flag) {
                             //reset_values(); 
                            print_labels();
                            first_locus_flag=0;
                        }

                           str_ptr=LOCUS_str;
                           BEGIN DEFINITION;
                        }

<DEFINITION>"DEFINITION"         {
                                      if( PRINT )  printf("%s\t",yytext);
                                     str_ptr=DEFINITION_str;
                                      BEGIN ACCESSION;
                                    }

<ACCESSION>"ACCESSION"         {    
                                      if( PRINT )  printf("%s\t",yytext);
                                      str_ptr=ACCESSION_str;
                                      BEGIN VERSION;
                                    }

<VERSION>"VERSION"         {       
                                if( PRINT )  printf("%s\t",yytext);
                                str_ptr=VERSION_str;
                                      BEGIN KEYWORDS;
                            }

<KEYWORDS>"KEYWORDS"         {
                                if( PRINT )  printf("%s\t",yytext);
                                    str_ptr=KEYWORDS_str;
                                      BEGIN SOURCE;
                                    }

<SOURCE>"SOURCE"         {   
                                if( PRINT )  printf("%s\t",yytext);
                                str_ptr=SOURCE_str;
                                      BEGIN ORGANISM;
                                    }

<ORGANISM>"ORGANISM"         { 
                                if( PRINT )  printf("%s\t",yytext);
                                str_ptr=ORGANISM_str;
                                BEGIN REFERENCE;
                             }

<REFERENCE>"REFERENCE"         {
                                if( PRINT )  printf("%s\t",yytext);
                                     str_ptr=REFERENCE_str;
                                      BEGIN AUTHORS;
                                    }

<AUTHORS>"AUTHORS"         {
                                if( PRINT )  printf("%s\t",yytext);
                                str_ptr=AUTHORS_str;
                                      BEGIN CONSRTM;
                                    }

<CONSRTM>"CONSRTM"         {      
                                if( PRINT )  printf("%s\t",yytext);
                                      str_ptr=CONSRTM_str;
                                      BEGIN TITLE;
                                    }

<TITLE>"TITLE"         {        
                                if( PRINT )  printf("%s\t",yytext);
                                   str_ptr=TITLE_str;
                                      BEGIN JOURNAL;
                                    }

<JOURNAL>"JOURNAL"         {    
                                if( PRINT )  printf("%s\t",yytext);
                                     str_ptr=JOURNAL_str;
                                      BEGIN PUBMED;
                                    }

<PUBMED>"PUBMED"         {          
                                if( PRINT )  printf("%s\t",yytext);
                                      str_ptr=PUBMED_str;
                                      BEGIN COMMENT;
                                    }

<COMMENT>"COMMENT"         {         
                                if( PRINT )  printf("%s\t",yytext);
                                      str_ptr=COMMENT_str;
                                      BEGIN FEATURES;
                                    }

<FEATURES>"FEATURES"         {   
                                if( PRINT )  printf("%s\t",yytext);
                                      str_ptr=FEATURES_str;
                                      BEGIN L_source;
                                    }

<L_source>"source"         {        
                                if( PRINT )  printf("%s\t",yytext);
                                      str_ptr=L_source_str;
                                      BEGIN F_strain;
                                    }

<F_strain>"\/strain="         {      
                                if( PRINT )  printf("%s\t",yytext);
                                      str_ptr=F_strain_str;
                                      BEGIN F_organism;
                                    }

<F_organism>"\/organism="         {
                                if( PRINT )  printf("%s\t",yytext);
                                      str_ptr=F_organism_str;
                                      BEGIN F_chromosome;
                                    }

<F_chromosome>"\/chromosome="         { 
                                if( PRINT )  printf("%s\t",yytext);
                                        str_ptr=F_chromosome_str;
                                      BEGIN gene;
                                    }

<gene>"gene "                        {
                                      if( PRINT )  printf("%s\t",yytext);
                                      partial_reset_values();
                                      gene_count++;
                                      //printf("Genecount = %d\n",gene_count);
                                      if( gene_count > MAX_GENES - 1  ) { 
                                           printf("Error: Too many genes. Increase the MAX_GENES definitioin\n");
                                           exit(0);
                                      }
                                      str_ptr=genes[gene_count].gene_str; 

                                      BEGIN gene_F_locus_tag;
                                    }

<gene_F_locus_tag>"\/locus_tag="      {
                                if( PRINT )  printf("%s\t",yytext);
                                      str_ptr=genes[gene_count].locus_tag;
                                      BEGIN CDS;
                                    }

<CDS>"gene "                        { 
                                if( PRINT )  printf("%s\t",yytext);
                                      gene_count++;
                                      if( gene_count > MAX_GENES - 1  ) { 
                                           printf("Error: Too many genes. Increase the MAX_GENES definitioin\n");
                                           exit(0);
                                      }
                                      str_ptr=genes[gene_count].gene_str; 
                                      BEGIN gene_F_locus_tag;
                                    }

<CDS>"CDS"                          {
                                if( PRINT )  printf("%s\t",yytext);
                                     str_ptr=genes[gene_count].CDS_str;
                                     /* BEGIN F_product; */
                                    }

<CDS>"\/product="                    {  
                                if( PRINT )  printf("%s\t",yytext);
                                        str_ptr=genes[gene_count].product;
                                   //   BEGIN F_codon_start;
                                    }

<CDS>"\/codon_start="                {  
                                if( PRINT )  printf("%s\t",yytext);
                                        str_ptr=genes[gene_count].codon_start;
                                  //    BEGIN F_locus_tag;
                                    }


<CDS>"\/locus_tag="         {           
                                if( PRINT )  printf("%s\t",yytext);
                                         str_ptr=F_locus_tag_str;
                       //               BEGIN F_transl_table;
                                    }


<CDS>"\/transl_table="    {             
                                if( PRINT )  printf("%s\t",yytext);
                                         str_ptr=genes[gene_count].transl_table;
                        //              BEGIN F_translation;
                                    }


<CDS>"\/note="                    {     
                                if( PRINT )  printf("%s\t",yytext);
                                        if( strlen(genes[gene_count].note) == 0 ) {
                                           str_ptr=genes[gene_count].note;
                                        }
                                        else {
                                            str_ptr = genes[gene_count].note + strlen(genes[gene_count].note);
                                        }
                                        
                         //             BEGIN F_db_xref;
                                    }

<CDS>"\/db_xref="         {             
                                if( PRINT )  printf("%s\t",yytext);
                                        if( strlen(genes[gene_count].db_xref) == 0 ) {
                                           str_ptr=genes[gene_count].db_xref;
                                        }
                                        else {
                                            str_ptr = genes[gene_count].db_xref + \
                                                       strlen(genes[gene_count].db_xref);
                                        }

                          //            BEGIN F_EC_number;
                                    }

<CDS>"\/EC_number="         {     
                                if( PRINT )  printf("%s\t",yytext);
                                   str_ptr=genes[gene_count].EC_number;
                           //       BEGIN F_translation;
                                }

<CDS>"\/translation="         {     
                                if( PRINT )  printf("%s\t",yytext);
                                   str_ptr=genes[gene_count].translation;
                       //           BEGIN ORIGIN;
                                }


<CDS>"ORIGIN"                     {   
                                if( PRINT )  printf("%s\t",yytext);
                                      str_ptr=ORIGIN_str;
                                      BEGIN NULL_STATE;
                                    }



<gene>"ORIGIN"                   { 
                                if( PRINT )  printf("%s\t",yytext);
                                  str_ptr=ORIGIN_str;
                            
                                   BEGIN NULL_STATE;
                          }



.       { 
            strncpy(str_ptr,yytext,strlen(yytext));
            str_ptr += strlen(yytext);
            *(str_ptr+1) ='\0';

         }
\n     {  }



%%
/*** C Code section ***/

int main(int argc, char **argv) 
{
     char *input, *output;
     int i;

     for(i = 1; i < argc && argc > 2 ; i++)
     {
        if( strncmp(argv[i], "-i", strlen("-i")) == 0 )
        {
            input = argv[++i];
        }
        else if( strncmp(argv[i], "-o", strlen("-o")) == 0 )
        {
            output     = argv[++i];
        }
       else {
            perror("Usage: ./genebank_parser -i inputfile -o outputfile" );
       }
     }

     printf("Gene count=%d\n",gene_count);
     FILE *infile = fopen(input, "r");
     if (!infile) {
        printf("Cannot  open  the annotation file %s\n",input);
        return(-1);
     }

     FILE *outfile = fopen(output, "w");
     if (!outfile) {
        printf("Cannot  open  the output file %s\n",input);
        return(-1);
    }


     yyin = infile;
     yyout= outfile;
     yylex();
     printf("Gene count=%d\n",gene_count);
     print_and_reset();
     return(0);
}


 /*    geneptr->codon_start[0]='\0';
     geneptr->locus_tag[0]='\0';
     geneptr->CDS[0]='\0';
     geneptr->product[0]='\0';
     geneptr->transl_table[0]='\0';
     geneptr->note[0]='\0';
     geneptr->db_xref[0]='\0';
     geneptr->EC_number[0]='\0';
     geneptr->translation[0]='\0';
*/


void partial_reset_values() {
     gene_str[0]='\0';
     F_codon_start_str[0]='\0';
     F_locus_tag_str[0]='\0';
     CDS_str[0]='\0';
     F_product_str[0]='\0';
     F_codon_str[0]='\0';
     F_transl_table_str[0]='\0';
     F_note_str[0]='\0';
     F_db_xref_str[0]='\0';
     F_EC_number_str[0]='\0';
     F_translation_str[0]='\0';
}

void reset_values() {
     LOCUS_str[0]='\0';
     DEFINITION_str[0]='\0';
     ACCESSION_str[0]='\0';
     VERSION_str[0]='\0';
     KEYWORDS_str[0]='\0';
     SOURCE_str[0]='\0';
     ORGANISM_str[0]='\0';
     REFERENCE_str[0]='\0';
     AUTHORS_str[0]='\0';
     CONSRTM_str[0]='\0';
     TITLE_str[0]='\0';
     JOURNAL_str[0]='\0';
     PUBMED_str[0]='\0';
     COMMENT_str[0]='\0';
     FEATURES_str[0]='\0';
     L_source_str[0]='\0';
     F_strain_str[0]='\0';
     F_organism_str[0]='\0';
     F_chromosome_str[0]='\0';

     gene_count=-1;
     int i;
     for(i=0; i< 10; i++) {
         genes[i].gene_str[0]='\0';
         genes[i].codon_start[0]='\0';
         genes[i].locus_tag[0]='\0';
         genes[i].CDS_str[0]='\0';
         genes[i].product[0]='\0';
         genes[i].transl_table[0]='\0';
         genes[i].note[0]='\0';
         genes[i].db_xref[0]='\0';
         genes[i].EC_number[0]='\0';
         genes[i].translation[0]='\0';
     }

     ORIGIN_str[0]='\0';
  
}

void print_value(const char *label, char *value) {
   fprintf(yyout, "<%s> : %s\n",label, value);
}

void print_value_delim(char *value, char *delim) {
   fprintf(yyout, "%s%s",value, delim);
}



