PROGS = mincmorph
OBJS = $(PROGS:=.o) kernel_io.o kernel_ops.o

CC=cc

OPTIONS = -O3 -w2
INCLUDES = -I/usr/local/include
CFLAGS = $(OPTIONS) $(INCLUDES)

LDINCLUDES = -L/usr/local/lib32
LDLIBS = -lvolume_io -lminc -lnetcdf -lm
LDOPTS = $(LDINCLUDES) $(LDLIBS)


all: $(PROGS) 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGS): $(OBJS)
	$(CC) $(OBJS) -o $@ $(OPTIONS) $(LDOPTS)

clean:
	rm -f *.o *~ $(PROGS)
