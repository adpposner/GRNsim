TARGETSIM = mainsim
TARGETNETGEN = mainnetgen

CC = gcc
export CC
#CFLAGS = -xHost -std=c99 -I/usr/include -ipo
CFLAGS = -Wunused -std=c99  
#LFLAGS =  -L/usr/local/lib -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lcrypto -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lxml2 
LFLAGS =  lib/mkl_small.so  -lm -ldl -lxml2 

#SRCS_NAMES = event minHeap reactionrates reaction genesListIO/reactantref genesListIO/basicgeneticelement genesListIO/globals genesListIO/cjson randomnumbers iofuncs parameters

SRCS_FILES = $(patsubst %,%.c,$(SRCS_NAMES))

OBJS_DIRS = sim models util iofuncs
OBJS_FILES = $(wildcard build/*/*.o)

export CFLAGS

all: CFLAGS += -O2 -DUNIFORM_TIME -DMESSENGERDECAYFORCESMICRODECAY -DPRINT_RNA
all: $(TARGETSIM)
all: $(TARGETNETGEN)
CMDLINEARGS =

recyc: CFLAGS += -DMESSENGERDECAYFORCESMICROUNBINDING -O2 -DUNIFORM_TIME -DPRINT_RNA
recyc: $(TARGETSIM)
recyc: $(TARGETNETGEN)

debug: CFLAGS += -g -O0 -DFEW_CYCLES -DDEBUG -DTESTASSERTIONS -DWITHASSERTIONS -DUNIFORM_TIME 
debug: $(TARGETSIM)
debug: $(TARGETNETGEN)

#Basic optimization flags
opt: CFLAGS += -O3 -xHost -DUNIFORM_TIME
opt: $(TARGETSIM)
opt: $(TARGETNETGEN)

#For Intel Advisor usage
advise: CFLAGS += -g -02 -vec -simd -DUNIFORM_TIME
advise: $(TARGETSIM)
advise: $(TARGETNETGEN)

#For profiling loops/functions
profile: CFLAGS += -profile-loops=all -O2 -DUNIFORM_TIME
profile: $(TARGETSIM)
profile: $(TARGETNETGEN)

o2only: CFLAGS += -DFEW_CYCLES -O2 -DUNIFORM_TIME
o2only: $(TARGET)

o3only: CFLAGS += -DFEW_CYCLES -O3 -DUNIFORM_TIME
o3only: $(TARGET)


subdirs:
	$(MAKE) -C iofuncs
	$(MAKE) -C util
	$(MAKE) -C models
	$(MAKE) -C sim

$(TARGETSIM): subdirs
	$(CC) $(CFLAGS) -o $(TARGETSIM) mainsim.c $(OBJS_FILES) $(LFLAGS)


$(TARGETNETGEN): subdirs
	$(CC) $(CFLAGS) -o $(TARGETNETGEN) mainnetgen.c $(OBJS_FILES) $(LFLAGS)



clean:
	rm $(OBJS_FILES) mainsim mainnetgen

cleanobjfiles:
	rm $(OBJS_FILES)

cleanoutput:
	rm -rf output/net*
