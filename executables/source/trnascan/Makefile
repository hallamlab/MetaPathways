# The location of the expat directory
CC=gcc  
LEX=lex  
LEXFLAGS=-lfl

all:  trnascan-1.4
	rm -rf trnascan.o

trnascan-1.4: trnascan.o
	$(CC)  -o trnascan-1.4 trnascan.c

clean:	
	rm -rf trnascan-1.4
#rm outline
#rm line
