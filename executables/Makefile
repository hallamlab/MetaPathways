# The location of the expat directory
CC=gcc  
LEX=lex  
LEXFLAGS=-lfl

all:  genbank2sequin lexCfile  trnascan-1.4
	$(CC) $(LEXFLAGS)  lex.yy.c  -o genbank_parser
	rm -rf lex.yy.c
	rm -rf trnascan.o

genbank2sequin:   genbank_to_sequin.l
	$(LEX)  $< 
	$(CC) $(LEXFLAGS)  lex.yy.c  -o genbank2sequin
	rm -rf lex.yy.c

lexCfile:   parse_genebank_to_table.l
	$(LEX)  parse_genebank_to_table.l 

trnascan-1.4: trnascan.o
	$(CC)  -o trnascan-1.4 trnascan.c

clean:	
	rm -rf genbank_parser trnascan-1.4 genbank2sequin
	rm -rf lex.yy.c
#rm outline
#rm line
