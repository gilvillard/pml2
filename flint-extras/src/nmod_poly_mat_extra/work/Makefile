INC_DIR = ../..
SYSLIB_DIR = /usr/local/lib
OBJS = $(wildcard $(LIB_DIR)/*.o)
SRC = $(wildcard test*c)
EXEC = $(patsubst %.c, %, $(SRC))
CHK = $(patsubst %.c, %.chk, $(patsubst test_%, %, $(SRC)))
COMPILE = $(patsubst %.c, %.exe, $(patsubst test_%, %, $(SRC)))

CC = gcc
CFLAGS = -Wall -O3 -march=native -I$(INC_DIR) -I/usr/local/include -I/opt/homebrew/include -I/Users/gvillard/include  -L/Users/gvillard/lib -L$(SYSLIB_DIR) -L/opt/homebrew/lib 
LIBS = -lflint -lmpfr -lpml -lgmp -lm -lpthread

BFLAGS = -Wall -O3 -march=native -I$(INC_DIR) -I/usr/local/include -I/opt/homebrew/include -I/Users/gvillard/include  -L/Users/gvillard/lib -L$(SYSLIB_DIR) -L/opt/homebrew/lib -lflint -lmpfr -lpml -lgmp -lm -lpthread

.PHONY: clean

clean:
	rm -f $(EXEC) $(CHK)

all:	$(COMPILE)

run:	$(CHK)

%.exe:
	if test -e test_$*.c; then $(CC) -o test_$* test_$*.c $(OBJS) $(CFLAGS) $(LIBS); fi

bench: bench.c 
	$(CC) -o bench bench.c $(BFLAGS)  

run: run.c 
	$(CC) -o run run.c $(BFLAGS)  

%.chk:
	if test -e test_$*; then ./test_$* > $*.chk; fi
