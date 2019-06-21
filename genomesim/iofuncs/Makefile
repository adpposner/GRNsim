CC = icc
BUILDDIR = $(patsubst %,%/iofuncs,$(BUILDDIR_BASE))


LOCALCFLAGS = -I/usr/local/include/graphviz
NEWLFLAGS =  -L/usr/local/lib -lgvc -lcgraph -lcdt

SRCS_NAMES =  json_interface iofuncs cjson

SRCS_FILES = $(patsubst %,%.c,$(SRCS_NAMES))

OBJS_FILES = $(patsubst %,$(BUILDDIR)/%.o,$(SRCS_NAMES))

CFLAGS += $(LOCALCFLAGS)

#all: appendlocalflags
all: makedirec $(OBJS_FILES)
	
makedirec:
	$(MKDIR_P) $(BUILDDIR)
	
$(BUILDDIR)/%.o: %.c
	icc $(CFLAGS) -c $^ -o $@ $(NEWLFLAGS)

#debug:	appendlocalflags
debug: $(OBJS_FILES)

#appendlocalflags: CFLAGS += LOCALCFLAGS
#	@echo appendlocalflags

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
	rm $(OBJS_FILES)

cleanexponentchange: ;