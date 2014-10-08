# Compute commutators from Trotter error term for a given n, where n
# is the number of orbitals. Symbolically stores both hamiltonian terms
# and coefficients that may be plugged in later.

CC = g++ -std=c++11
CFLAGS = -Wall -g
LDFLAGS = -lm

all : compute-commutators

compute-commutators-util.o : compute-commutators-util.cc \
	compute-commutators-util.h
	${CC} ${CFLAGS} -c compute-commutators-util.cc

compute-commutators.o : compute-commutators.cc compute-commutators.h \
	compute-commutators-util.h
	${CC} ${CFLAGS} -c compute-commutators.cc compute-commutators-util.cc

compute-commutators_main.o : compute-commutators_main.cc compute-commutators.h 
	${CC} ${CFLAGS} -c compute-commutators_main.cc

compute-commutators : compute-commutators-util.o compute-commutators.o \
	compute-commutators_main.o
	${CC} ${CFLAGS} compute-commutators-util.o compute-commutators.o \
	compute-commutators_main.o -o compute-commutators

clean :
	rm -rf *.o all
