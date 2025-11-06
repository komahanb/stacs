include ${HOME}/git/tacs/Makefile.in
include ${HOME}/git/tacs/TACS_Common.mk

PSPACE_INCLUDE = -I${HOME}/git/pspace/src/include
PSPACE_LIB = ${HOME}/git/pspace/cpp/libpspace.a

PSPACE_TACS_INCLUDE = -I${HOME}/git/pspace/examples/tacs/
PSPACE_TACS_LIB = ${HOME}/git/pspace/examples/tacs/libtacsuq.a

# This is the one rule that is used to compile all the
# source code in TACS
%.o: %.cpp
	${CXX} ${TACS_CC_FLAGS} ${PSPACE_INCLUDE} ${PSPACE_TACS_INCLUDE} -I../../src/include -I/usr/include/python3.10 -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.cpp successfully ---"
	@echo

OBJS = deterministic.o sampling.o projection.o

default: ${OBJS}
	${CXX} -o deterministic deterministic.o ${TACS_LD_FLAGS} ${PSPACE_LIB} ${PSPACE_TACS_LIB} -L../../lib -lstacs
	${CXX} -o sampling sampling.o ${TACS_LD_FLAGS} ${PSPACE_LIB} ${PSPACE_TACS_LIB} -L../../lib -lstacs
	${CXX} -o projection projection.o ${TACS_LD_FLAGS} ${PSPACE_LIB} ${PSPACE_TACS_LIB} -L../../lib -lstacs

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF=-DTACS_USE_COMPLEX -DUSE_COMPLEX
complex: default

complex_debug: TACS_DEF=-DTACS_USE_COMPLEX -DUSE_COMPLEX
complex_debug: debug

clean:
	rm -f *.o deterministic sampling projection

test: default
	./deterministic
	./sampling
	./projection

test_complex: complex
	./deterministic
	./sampling
	./projection
