BUILDDIR = ../build/iofuncs

LOCALCFLAGS =
LFLAGS =
SRCS_NAMES =  json_interface iofuncs cjson

SRCS_FILES = $(patsubst %,%.c,$(SRCS_NAMES))

OBJS_FILES = $(patsubst %,$(BUILDDIR)/%.o,$(SRCS_NAMES))

CFLAGS += $(LOCALCFLAGS)


all: $(OBJS_FILES)
	

$(BUILDDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(LFLAGS)

debug: $(OBJS_FILES)

opt:
	$(CC) -O3 -xHost $(CFLAGS) -o $(TARGET) main.c $(SRCS_FILES) $(LFLAGS)


clean:
	rm $(OBJS_FILES)


#Associated tests - all stored in geneslistTests.c originally
# testJSON:
# 	$(CC) -g -O0 -DTESTREACTANTREF $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS)

# testRxns:	
# 	$(CC) -g -O0 -DTESTREACTIONGENERATION $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS) 

# testConnections:
# 	$(CC) -g -O0 -DTESTPTRCONVERSION $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS)

# testInit:
# 	$(CC) -g -O0 -DTESTASSERTIONS -DCONVERSIONANDINIT $(CFLAGS) -o $(TARGET) geneslistTests.c $(SRCS_FILES) $(NEWLFLAGS)
