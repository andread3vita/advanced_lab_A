CXX=g++
CFLAGS=-Wall -fPIC -shared
SRC=src
LIB=lib
SRCS=$(wildcard $(SRC)/*.cc)
LIBS=$(patsubst $(SRC)/%.cc, $(LIB)/lib%.so, $(SRCS))


all:$(LIBS)

$(LIB)/lib%.so: $(SRC)/%.cc
	$(CXX) $(CFLAGS) -o $@ $^ 

