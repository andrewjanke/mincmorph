PROGS = mincmorph
HEADERS = kernel_io.h kernel_ops.h
OBJS = $(PROGS:=.o) $(HEADERS:.h=.o)

CC=cc

OPTIONS = -g3 -O3 -fullwarn -w2
INCLUDES = -I/usr/local/include
CFLAGS = $(OPTIONS) $(INCLUDES)

LDINCLUDES = -L/usr/local/lib32
LDLIBS = -lvolume_io -lminc -lnetcdf -lm
LDOPTS = $(LDINCLUDES) $(LDLIBS)


all: $(PROGS)

.c.o: $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGS): $(OBJS)
	$(CC) $(OBJS) -o $@ $(OPTIONS) $(LDOPTS)

clean:
	rm -f *.o *~ $(PROGS)
