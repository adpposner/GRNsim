CC = icc
MYBUILDDIR = netgen
#BUILDDIR = $(patsubst %,$(BUILDDIR_BASE)/%,$(MYBUILDDIR))
BUILDDIR = $(patsubst %,%/netgen,$(BUILDDIR_BASE))

LOCALCFLAGS = -I/usr/include/libxml2
LFLAGS = -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -lxml2

CFLAGS += $(LOCALCFLAGS)
SRCS_NAMES = networkdistributions

SRCS_FILES = $(patsubst %,%.c,$(SRCS_NAMES))

OBJS_FILES = $(patsubst %,$(BUILDDIR)/%.o,$(SRCS_NAMES))

all: makedirec $(OBJS_FILES)
#all: $(OBJS_FILES)
	
makedirec:
	$(MKDIR_P) $(BUILDDIR)

$(BUILDDIR)/%.o: %.c
	icc $(CFLAGS) -c $^ -o $@ $(NEWLFLAGS)


debug: $(OBJS_FILES)

verbose: CFLAGS += -DDEBUG
verbose: debug

opt:
	icc -O3 -xHost $(CFLAGS) -o $(TARGET) main.c $(SRCS_FILES) $(NEWLFLAGS)

clean:
	rm $(BUILDDIR)/*.o

cleanexponentchange: ;
