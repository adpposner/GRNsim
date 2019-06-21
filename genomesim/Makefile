TARGETSIM = mainsim
TARGETNETGEN = mainnetgen
BUILDDIR_ROOT = build2
BUILDDIR_BASE = $(addprefix ../,$(BUILDDIR_ROOT))
export BUILDDIR_BASE
MKDIR_P = mkdir -p
export MKDIR_P
CC = icc

CFLAGS = -xHost -std=c99 -I/usr/include/graphviz -I/usr/include -I/opt/local/include 
LFLAGS = -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

NEWLFLAGS =  -L/usr/local/lib -L${MKLROOT}/lib/intel64 -Wl,-rpath,${MKLROOT}/lib/intel64 -lcrypto -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lxml2 

FORCES_DECAY = -DMESSENGERDECAYFORCESMICRODECAY
CMDLINEARGS = 
CFLAGS += $(CMDLINEARGS)
#SRCS_NAMES = event minHeap reactionrates reaction graph_interface genesListIO/reactantref genesListIO/basicgeneticelement genesListIO/globals genesListIO/cjson randomnumbers iofuncs parameters

SRCS_FILES = $(patsubst %,%.c,$(SRCS_NAMES))

OBJS_DIRS = sim netgen models iofuncs util
OBJS_FILES = $(wildcard $(BUILDDIR_ROOT)/*/*.o)

export CFLAGS
all: makedirec
all: CFLAGS += -O2 $(FORCES_DECAY) -DUNIFORM_TIME -ipo
all: $(TARGETSIM)
#all: cleanofiles
all: $(TARGETNETGEN)

.PHONY: subdirs $(OBJS_DIRS) cleanexponentchange

makedirec:
	$(MKDIR_P) $(BUILDDIR_ROOT)

subdirs: $(OBJS_DIRS)

$(OBJS_DIRS):
	$(MAKE) -C $@

$(TARGETSIM): subdirs
	$(CC) $(CFLAGS) -o $(TARGETSIM) simmain.c $(OBJS_FILES) $(NEWLFLAGS)



$(TARGETNETGEN): subdirs
	$(CC) $(CFLAGS) -o $(TARGETNETGEN) networkgenmain.c $(OBJS_FILES) $(NEWLFLAGS)


#MESSENGERDECAYFORCESUNBINDING
debug: CFLAGS += -g -O0 -DFEW_CYCLES -DDEBUG $(FORCES_DECAY) -DTESTASSERTIONS -DWITHASSERTIONS -DUNIFORM_TIME 
debug: $(TARGETSIM)
debug: $(TARGETNETGEN)

verbose: CFLAGS += -DDEBUG
verbose: debug

advise: CFLAGS += -g -02 -vec -simd $(FORCES_DECAY) -DUNIFORM_TIME -ipo
advise: $(TARGETSIM)
advise: $(TARGETNETGEN)

pgoA: CFLAGS += -O3 -xHost $(FORCES_DECAY) -DUNIFORM_TIME  -DPRINT_RNA -prof-gen -prof-dir=/home/russell/Desktop/cqmprojs/genomesim/pgo
pgoA: $(TARGETSIM)

pgoB: CFLAGS += -O3 -xHost $(FORCES_DECAY) -DUNIFORM_TIME -DPRINT_RNA -prof-use -prof-dir=/home/russell/Desktop/cqmprojs/genomesim/pgo -ipo
pgoB: $(TARGETSIM)


opt: CFLAGS += -O3 -xHost $(FORCES_DECAY) -DUNIFORM_TIME -ipo
opt: $(TARGETSIM)
opt: $(TARGETNETGEN)

o2only: CFLAGS += -DFEW_CYCLES -O2 $(FORCES_DECAY) -DUNIFORM_TIME -ipo
o2only: $(TARGET)

o3only: CFLAGS += -DFEW_CYCLES -O3 $(FORCES_DECAY) -DUNIFORM_TIME -ipo
o3only: $(TARGET)

profile: CFLAGS += -profile-loops=all -O2 $(FORCES_DECAY) -DUNIFORM_TIME -ipo
profile: $(TARGETSIM)
profile: $(TARGETNETGEN)



cleanexponentchange:
	rm build/sim/reactionrates.o build/models/basicgeneticelement.o build/util/globals.o

clean:
	rm $(OBJS_FILES) mainsim mainnetgen

cleanofiles:
	rm $(OBJS_FILES)

cleanoutput:
	rm -rf output/net-*
