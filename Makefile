
# build options
OPTS = -DLOOP_CNT=1

# GCC
CC 	= gcc
CFLAGS1 = -DGCC -Wall
CFLAGS2 = -DGCC -Wall -funroll-loops -funroll-all-loops -O2
CFLAGS3 = -DGCC -Wall -funroll-loops -funroll-all-loops -O3

# Intel ICX
#CC      = icx
#CFLAGS1 = -DICX -Wall
#CFLAGS2 = -DICX -Wall -O2
#CFLAGS3 = -DICX -Wall -O3

#all: tc tcO2 tcO3
all: tcO2

tc: tc.c
	${CC} ${CFLAGS} ${CFLAGS1} ${OPTS} -o tc tc.c -lm

tcO2: tc.c
	${CC} ${CFLAGS} ${CFLAGS2} ${OPTS} -o tcO2 tc.c -lm

tcO3: tc.c
	${CC} ${CFLAGS} ${CFLAGS3} ${OPTS} -o tcO3 tc.c -lm

testO2: tcO2 roadNet-PA.mtx
	./tcO2 -f roadNet-PA.mtx -x

clean: 
	rm -f core *~ \
	tc tcO2 tcO3
