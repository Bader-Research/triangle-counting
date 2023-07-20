
# build options
#OPTS = -DLOOP_CNT=1

# Turn off/on OpenMP parallel code
#PARALLEL = 
PARALLEL = -DPARALLEL

WARN =
#WARN = -Wall

# GCC
CC 	= gcc
CFLAGS1 = -DGCC $(PARALLEL) $(WARN)
CFLAGS2 = -DGCC $(PARALLEL) $(WARN) -funroll-loops -funroll-all-loops -O2
CFLAGS3 = -DGCC $(PARALLEL) $(WARN) -funroll-loops -funroll-all-loops -O3

# Intel ICX
#CC      = icx
#CFLAGS1 = -DICX -DPARALLEL $(WARN)
#CFLAGS2 = -DICX -DPARALLEL $(WARN) -O2
#CFLAGS3 = -DICX -DPARALLEL $(WARN) -O3

OBJS = queue.o graph.o tc.o tc_parallel.o

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
