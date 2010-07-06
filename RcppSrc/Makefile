include ${R_HOME}/etc${R_ARCH}/Makeconf
SRC = Rcpp.cpp
OBJ = Rcpp.o
LIB = libRcpp.a

.PHONY: all

all: distclean $(LIB) clean

$(LIB): $(OBJ)
	$(AR) crs $(LIB) $(OBJ)

$(OBJ): $(SRC)

clean:
	$(RM) $(OBJ)

distclean:
	$(RM) $(OBJ) $(LIB)
