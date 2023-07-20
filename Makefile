
# build options
#OPTS = -DLOOP_CNT=1

# Turn off/on OpenMP parallel code
#PARALLEL = 
PARALLEL = -DPARALLEL

# GCC
CC 	= gcc
CFLAGS1 = -DGCC $(PARALLEL) -Wall
CFLAGS2 = -DGCC $(PARALLEL) -Wall -funroll-loops -funroll-all-loops -O2
CFLAGS3 = -DGCC $(PARALLEL) -Wall -funroll-loops -funroll-all-loops -O3

# Intel ICX
#CC      = icx
#CFLAGS1 = -DICX -DPARALLEL -Wall
#CFLAGS2 = -DICX -DPARALLEL -Wall -O2
#CFLAGS3 = -DICX -DPARALLEL -Wall -O3

OBJS = tc.o tc_parallel.o

#all: tc tcO2 tcO3
all: tcO2

%.o: %.c
	${CC} ${CFLAGS} ${CFLAGS2} ${OPTS} -c $< -o $@

tc: $(OBJS)
	${CC} ${CFLAGS} ${CFLAGS1} ${OPTS} -o $@ $(OBJS) -lm

tcO2: $(OBJS)
	${CC} ${CFLAGS} ${CFLAGS2} ${OPTS} -o $@ $(OBJS) -lm

tcO3: $(OBJS)
	${CC} ${CFLAGS} ${CFLAGS3} ${OPTS} -o $@ $(OBJS) -lm

clean: 
	rm -f core *~ $(OBJS) \
	tc tcO2 tcO3
