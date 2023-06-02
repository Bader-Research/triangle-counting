
# GCC
CC 	= gcc
CFLAGS1 = -DGCC -Wall
CFLAGS2 = -DGCC -Wall -funroll-loops -funroll-all-loops -O2
CFLAGS3 = -DGCC -Wall -funroll-loops -funroll-all-loops -O3

#all: tc tcO2 tcO3
all: tcO2

tc: tc.c
	${CC} ${CFLAGS1} -o tc tc.c -lm

tcO2: tc.c
	${CC} ${CFLAGS2} -o tcO2 tc.c -lm

tcO3: tc.c
	${CC} ${CFLAGS3} -o tcO3 tc.c -lm

clean: 
	rm -f core *~ \
	tc tcO2 tcO3
