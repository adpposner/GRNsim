CC = icc
BUILDDIR = $(patsubst %,%/models,$(BUILDDIR_BASE))


CTESTFLAGS = -I/usr/include/libxml2 -I/usr/local/include/graphviz
#LFLAGS = -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

NEWLFLAGS =  -L/usr/local/lib -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lgvc -lcgraph -lcdt -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lxml2 

SRCS_NAMES =  minheap reaction geneslist basicgeneticelement 
#minheap reaction geneslist basicgeneticelement


SRCS_FILES = $(patsubst %,%.c,$(SRCS_NAMES))

OBJS_FILES = $(patsubst %,$(BUILDDIR)/%.o,$(SRCS_NAMES))

LOCALCFLAGS =
CFLAGS += $(LOCALCFLAGS)

all: makedirec $(OBJS_FILES)

	
makedirec:
	$(MKDIR_P) $(BUILDDIR)
#debug: CFLAGS += -DTESTASSERTIONS -O0 -g
debug:	$(OBJS_FILES)	


$(BUILDDIR)/%.o: %.c
	icc $(CFLAGS) -c $^ -o $@  $(NEWLFLAGS)

testmodels:
	icc -g -O0 -DTESTMODELS -DTESTASSERTIONS $(CTESTFLAGS) -o testme modelgenerator.c reaction.c ../iofuncs/json_interface.c geneslist.c basicgeneticelement.c ../iofuncs/cjson.c ../util/globals.c ../util/parameters.c ../util/randomnumbers.c $(NEWLFLAGS)

testRxns:	
	icc -g -O0 -DTESTREACTIONGENERATION $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS) 

testConnections:
	icc -g -O0 -DTESTPTRCONVERSION $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS)

testInit:
	icc -g -O0 -DTESTASSERTIONS -DCONVERSIONANDINIT $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS)

clean:
	rm ../build/models/*.o

cleanexponentchange:
	rm ../build/models/basicgeneticelement.o