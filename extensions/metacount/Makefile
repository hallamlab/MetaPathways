#CC = /home/sgeadmin/gcc_install/bin/g++ -std=c++17 -Wall # -g
CC = g++ -std=c++17 -Wall # -g
#CC = g++ -pg
#CC = g++ -Wall 
CCFLAGS=  -O4 -m64

CFLAGS=-Iparser/include -Lparser

PROG = metacount
SOURCES= utilities.cpp options.cpp  metacount.cpp core.cpp parser.cpp types.cpp

OBJECTS= $(SOURCES:.cpp=.o)
HEADERS= $(SOURCES:.cpp=.h)

LIB = parser/libbamparser.a
#LIB_PATH = /home/sgeadmin/gcc_install/lib64/

all: $(PROG)

%.o:%.cpp   $(SOURCES) types.h
	$(CC) $(CCFLAGS) $(CFLAGS) $< -c -o $@  

$(LIB):
	make -C parser

clean:
	make -C parser clean
	rm -rf $(OBJECTS) $(PROG)

$(PROG): $(LIB) $(OBJECTS) $(HEADERS) types.h 
	$(CC) $(CCFLAGS) $(CFLAGS)  $(OBJECTS) -o $(PROG) -lbamparser -lz -lstdc++fs -static

