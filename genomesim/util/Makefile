BUILDDIR = ../build/util


LOCALCFLAGS = -I/usr/include/libxml2 -I/opt/local/include

LFLAGS =  -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lcrypto -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lxml2 

SRCS_NAMES = modelloader modelgenerator parameters globals randomnumbers

SRCS_FILES = $(patsubst %,%.c,$(SRCS_NAMES))

OBJS_FILES = $(patsubst %,$(BUILDDIR)/%.o,$(SRCS_NAMES))

CFLAGS += $(LOCALCFLAGS)

all: $(OBJS_FILES)
	

$(BUILDDIR)/%.o: %.c
	$(CC) $(CFLAGS)  -c $^ -o $@ $(LFLAGS)

debug: $(OBJS_FILES)

opt:
	$(CC) -O3 -xHost $(CFLAGS) -o $(TARGET) main.c $(SRCS_FILES) $(LFLAGS)


clean:
	rm ../build/util/*.o
