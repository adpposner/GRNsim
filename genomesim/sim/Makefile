CC = icc
BUILDDIR = $(patsubst %,%/sim,$(BUILDDIR_BASE))


LOCALCFLAGS = 
LFLAGS = -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

CFLAGS += $(LOCALCFLAGS)
SRCS_NAMES = simulator event reactionrates

SRCS_FILES = $(patsubst %,%.c,$(SRCS_NAMES))

OBJS_FILES = $(patsubst %,$(BUILDDIR)/%.o,$(SRCS_NAMES))


all: makedirec $(OBJS_FILES)


makedirec:
	$(MKDIR_P) $(BUILDDIR)

$(BUILDDIR)/%.o: %.c 
	$(CC) $(CFLAGS) -c $^ -o $@ $(NEWLFLAGS)



debug: $(OBJS_FILES)

verbose: CFLAGS += -DDEBUG
verbose: debug

opt:
	icc -O3 -xHost $(CFLAGS) -o $(TARGET) main.c $(SRCS_FILES) $(NEWLFLAGS)

testJSON:
	icc -g -O0 -DTESTREACTANTREF $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS)

testRxns:	
	icc -g -O0 -DTESTREACTIONGENERATION $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS) 

testConnections:
	icc -g -O0 -DTESTPTRCONVERSION $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS)

testInit:
	icc -g -O0 -DTESTASSERTIONS -DCONVERSIONANDINIT $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS)

clean:
	rm $(wildcard $(BUILDDIR)/*.o)

cleanexponentchange:
	rm ../build/sim/reactionrates.o
