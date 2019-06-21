BUILDDIR = ../build/sim



LOCALCFLAGS = 
LFLAGS = -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

CFLAGS += $(LOCALCFLAGS)
SRCS_NAMES = simulator event reactionrates

SRCS_FILES = $(patsubst %,%.c,$(SRCS_NAMES))

OBJS_FILES = $(patsubst %,$(BUILDDIR)/%.o,$(SRCS_NAMES))


all: $(OBJS_FILES)
	

$(BUILDDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(NEWLFLAGS)


debug: $(OBJS_FILES)

verbose: CFLAGS += -DDEBUG
verbose: debug

opt:
	$(CC) -O3 -xHost $(CFLAGS) -o $(TARGET) main.c $(SRCS_FILES) $(NEWLFLAGS)

testJSON:
	$(CC) -g -O0 -DTESTREACTANTREF $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS)

testRxns:	
	$(CC) -g -O0 -DTESTREACTIONGENERATION $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS) 

testConnections:
	$(CC) -g -O0 -DTESTPTRCONVERSION $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS)

testInit:
	$(CC) -g -O0 -DTESTASSERTIONS -DCONVERSIONANDINIT $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS)

clean:
	rm ../build/sim/*.o
