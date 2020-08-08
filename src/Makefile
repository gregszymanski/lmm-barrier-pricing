.PHONY = clean, mrproper

DEBUG = no
CC = g++ 
BOOST = yes
EIGEN = yes

EXEC = barrier_option_libor



ifeq ($(DEBUG),yes)
    CFLAGS = -g -ansi -Wall -pedantic -Wextra -DDEBUG -std=gnu++14 -Wno-c99-extensions 
    LDFLAGS = -lm 
else
    CFLAGS = -std=gnu++14 -O2
    LDFLAGS = -lm 
endif

ifeq ($(EIGEN), yes)
    EIGENLIB = -I ./libs/Eigen
else
    EIGENLIB = 
endif

ifeq ($(BOOST), yes)
    BOOSTLIB = -I ./libs/boost_1_71_0
else
    BOOSTLIB = 
endif



FILE = $(wildcard *.cpp)

FILE_OBJ = $(FILE:.cpp=.o)




all: $(EXEC)
	@cp $(EXEC) ..
	@echo "Compilation du projet"

install: $(EXEC)
	@cp $(EXEC) ..
	@echo "installation du projet"


$(EXEC): $(FILE_OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.cpp 
	$(CC) $(BOOSTLIB) $(EIGENLIB) -o $@ -c $< $(CFLAGS)
	
clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)


