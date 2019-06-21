CC=icc
FORTC=ifort
TARGET=main


IFORTFLAGSLINUX=-I${MKLROOT}/include
LFORTFLAGSLINUX=-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

IFORTFLAGSMAC=-I${MKLROOT}/include
LFORTFLAGSMAC=-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

CFLAGS=
FORTFLAGS=

UNAME := $(shell uname)

C_SRC=main.c
F_SRC=htmod.f90
C_OBJ=main.o
F_OBJ=htmod.o

ifeq ($(UNAME),Linux)
IFORTFLAGS=$(IFORTFLAGSLINUX)
LFORTFLAGS=$(LFORTFLAGSLINUX)
else 
IFORTFLAGS=$(IFORTFLAGSMAC)
LFORTFLAGS=$(LFORTFLAGSMAC)
endif

FORTFLAGS=
CFLAGS=


all: $(TARGET)

debug: FORTFLAGS+=-g -traceback
debug: CFLAGS+=-g
debug: $(TARGET)

$(TARGET): CPART FPART
	$(FORTC) $(FORTFLAGS) -o $(TARGET) -nofor-main $(F_OBJ) $(C_OBJ) $(LFORTFLAGS)

mainfort: FPART
	$(FORTC) $(FORTFLAGS) -o mainfort mainfort.f90 $(F_OBJ) $(LFORTFLAGS)


CPART:
	$(CC) $(CFLAGS) -c $(C_SRC)

FPART:
	$(FORTC) $(FORTFLAGS) -c $(F_SRC) $(IFORTFLAGS) $(LFORTFLAGS)

.PHONY: clean
clean:
	rm *.o $(TARGET)
