CC = icc
BUILDDIR = $(patsubst %,%/util,$(BUILDDIR_BASE))


LOCALCFLAGS = -I/usr/include/libxml2
LFLAGS = -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm 

NEWLFLAGS =  -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lcrypto -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lxml2 

SRCS_NAMES = modelloader modelgenerator parameters globals randomnumbers

SRCS_FILES = $(patsubst %,%.c,$(SRCS_NAMES))

OBJS_FILES = $(patsubst %,$(BUILDDIR)/%.o,$(SRCS_NAMES))

CFLAGS += $(LOCALCFLAGS)

all: makedirec $(OBJS_FILES)


makedirec:
	$(MKDIR_P) $(BUILDDIR)

$(BUILDDIR)/%.o: %.c
	icc $(CFLAGS)  -c $^ -o $@ $(NEWLFLAGS)

debug: $(OBJS_FILES)

opt:
	icc -O3 -xHost $(CFLAGS) -o $(TARGET) main.c $(SRCS_FILES) $(NEWLFLAGS)


clean:
	rm ../build/util/*.o

cleanexponentchange:
	rm ../build/util/globals.o