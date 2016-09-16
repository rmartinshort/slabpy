ARCH = x86_64

BDIR = bin/
ODIR = objects/
CFLAGS = -v3 -Wall


all: dirs $(BDIR)/asciiint2bebin $(BDIR)/asciifloat2bebin $(BDIR)/sphere $(BDIR)/sumfloatbin

dirs:
	if [ ! -s $(ODIR) ];then\
		mkdir -p $(ODIR); \
	fi;\
	if [ ! -s $(BDIR) ];then\
		mkdir -p $(BDIR); \
	fi

$(BDIR)/asciiint2bebin: asciiint2bebin.c $(ODIR)/flip_byte_order.o
	$(CC) $(CFLAGS) asciiint2bebin.c  $(ODIR)/flip_byte_order.o -o $(BDIR)/asciiint2bebin $(LDFLAGS)

$(BDIR)/asciifloat2bebin: asciifloat2bebin.c $(ODIR)/flip_byte_order.o
	$(CC) $(CFLAGS) asciifloat2bebin.c  $(ODIR)/flip_byte_order.o -o $(BDIR)/asciifloat2bebin $(LDFLAGS)

$(BDIR)/sphere: sphere.c
	$(CC) $(CFLAGS) sphere.c -o $(BDIR)/sphere -lm

$(BDIR)/sumfloatbin: sumfloatbin.c
	$(CC) $(CFLAGS) sumfloatbin.c -o $(BDIR)/sumfloatbin -lm

$(ODIR)/%.o: %.c  $(HDR_FLS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $(ODIR)/$*.o

