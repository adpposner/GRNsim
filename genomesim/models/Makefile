BUILDDIR = ../build/models

LFLAGS =  -L/usr/local/lib -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lxml2 

SRCS_NAMES =  minheap reaction geneslist basicgeneticelement ../iofuncs/json_interface


SRCS_FILES = $(patsubst %,%.c,$(SRCS_NAMES))

OBJS_FILES = $(patsubst %,$(BUILDDIR)/%.o,$(SRCS_NAMES))

LOCALCFLAGS =
CFLAGS += $(LOCALCFLAGS)


all: $(OBJS_FILES)

debug:	$(OBJS_FILES)	

$(BUILDDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@  $(LFLAGS)


clean:
	rm ../build/models/*.o


#Associated tests invoked by #defines
# testmodels:
# 	$(CC) -g -O0 -DTESTMODELS -DTESTASSERTIONS $(CTESTFLAGS) -o testme modelgenerator.c reaction.c ../iofuncs/json_interface.c geneslist.c basicgeneticelement.c ../iofuncs/cjson.c ../util/globals.c ../util/parameters.c ../util/randomnumbers.c $(LFLAGS)

# testRxns:	
# 	$(CC) -g -O0 -DTESTREACTIONGENERATION $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(LFLAGS) 

# testConnections:
# 	$(CC) -g -O0 -DTESTPTRCONVERSION $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(LFLAGS)

# testInit:
# 	$(CC) -g -O0 -DTESTASSERTIONS -DCONVERSIONANDINIT $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(LFLAGS)
