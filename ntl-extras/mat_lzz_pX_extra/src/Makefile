
INC_DIR = ../../include
LIB_DIR = ../../lib

SYS_INC = /usr/local/include

SOURCES = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp, $(LIB_DIR)/%.o, $(SOURCES))

THRESHOLDS_HEADERS = 

CXX = g++
CXXFLAGS = -Wall -std=c++11 -O3 -march=native -I$(INC_DIR) -I$(SYS_INC) -g
LIBS = -lntl -lgmp -lm
.PHONY: clean

all:	cp $(OBJS)

cp:
	for file in $(THRESHOLDS_HEADERS) ; do \
		if [ ! -f ../$$file.h ]; then cp -f ../$$file.orig $(INC_DIR)/$$file.h ; fi ; \
	done
	cp -f ../*.h $(INC_DIR)

clean:
	rm -f $(OBJS)

reset:	
	for file in $(THRESHOLDS_HEADERS) ; do \
		rm -f ../$$file.h; \
	done


$(LIB_DIR)/%.o:	%.cpp
	$(CXX) -c $< $(CXXFLAGS) -o $@

