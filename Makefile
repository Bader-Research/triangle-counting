
# build options
# OPTS = -DLOOP_CNT=1

# Turn off/on OpenMP parallel code
#PARALLEL = 
#GCC
PARALLEL = -DPARALLEL -fopenmp
#Intel ICX
#PARALLEL = -DPARALLEL -qopenmp

WARN =
#WARN = -Wall

# GCC
CC     = gcc
CFLAGS = -DGCC $(PARALLEL) $(WARN) -funroll-loops -funroll-all-loops -O2

# Intel ICX
#CC     = icx
#CFLAGS = -DICX $(PARALLEL) $(WARN) -O2

SRCS = $(wildcard *.c)
OBJS = $(SRCS:.c=.o)

all: tc

%.o: %.c
	${CC} ${CFLAGS} ${OPTS} -c $< -o $@

tc: $(OBJS)
	${CC} ${CFLAGS} ${OPTS} -o $@ $(OBJS) -lm

clean: 
	rm -f core *~ $(OBJS) tc
