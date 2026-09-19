.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# === Directory layout ========================================================
objdir  = obj
srcdir  = src
testdir = tests
bindir  = bin

# === Command-line switches ===================================================
debug  ?= 0
openmp ?= 1

defaulttarget: usage

# === Compiler configuration ==================================================
# Inserted by `python config.py config/<machine>_<compiler>`. The fragment
# defines the compiler (FC), flag sets (FFLAGS_BASE, FFLAGS_OPENMP, MODFLAGS,
# CPPFLAGS_PP, DFLAGS_*), the netCDF wiring (INC_NC, LIB_NC), and any
# host-specific link extras (LFLAGS_EXTRA).

<COMPILER_CONFIGURATION>

# === Shared dependency wiring ================================================
# Pulled in after the compiler fragment so it can reference the variables above.
# Defines the fesm-utils / FFTW / SHTns paths and the final FFLAGS_FE / LFLAGS_FE
# / CPPFLAGS_FE flag sets, honouring the openmp= switch.
include config/common.mk

# === Git version string ======================================================
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags 2>/dev/null || echo "unknown")
CPPFLAGS_FE := $(CPPFLAGS_FE) -DVERSION=\"$(GIT_VERSION)\"

# === Debug / profile selection ===============================================
DFLAGS = $(DFLAGS_NODEBUG)
ifeq ($(debug),1)
	DFLAGS = $(DFLAGS_DEBUG)
endif
ifeq ($(debug),2)
	DFLAGS = $(DFLAGS_PROFILE)
endif

# === Assemble the flags actually used by the rules ===========================
FFLAGS   = $(FFLAGS_FE)
LFLAGS   = $(LFLAGS_FE)
CPPFLAGS = $(CPPFLAGS_FE)
ifeq ($(openmp),1)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPENMP)
endif

# === Source/rule lists and build targets =====================================
include config/Makefile_fastearth.mk
