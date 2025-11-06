include ${HOME}/git/tacs/Makefile.in
include ${HOME}/git/tacs/TACS_Common.mk

PSPACE_INCLUDE = -I${HOME}/git/pspace/src/include
PSPACE_LIB = ${HOME}/git/pspace/lib/libpspace.a

%.o: %.cpp
	${CXX} ${TACS_CC_FLAGS} ${PSPACE_INCLUDE} ${PSPACE_TACS_INCLUDE} -I../../src/include -I/usr/include/python3.10 -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.cpp successfully ---"
	@echo

COMMON_OBJS = common.o TACSPotentialEnergy.o TACSDisplacement.o

default: deterministic sampling projection

deterministic: ${COMMON_OBJS} deterministic.o
	${CXX} -o $@ ${COMMON_OBJS} deterministic.o ${TACS_LD_FLAGS} -L../../lib -lstacs -Wl,-rpath,'$$ORIGIN/../../lib' ${PSPACE_LIB}

sampling: ${COMMON_OBJS} sampling.o
	${CXX} -o $@ ${COMMON_OBJS} sampling.o ${TACS_LD_FLAGS} -L../../lib -lstacs -Wl,-rpath,'$$ORIGIN/../../lib' ${PSPACE_LIB}

projection: ${COMMON_OBJS} projection.o
	${CXX} -o $@ ${COMMON_OBJS} projection.o ${TACS_LD_FLAGS} -L../../lib -lstacs -Wl,-rpath,'$$ORIGIN/../../lib' ${PSPACE_LIB}

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
