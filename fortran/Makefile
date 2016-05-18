SHELL = /bin/sh
DEBUG = false
FC = gfortran
CC = gcc
FCFLAGS = -g -O2 -Wall

#Adapt to host system
ifeq (,$(findstring Windows,$(OS)))
  RM = rm -f
else
  RM = del
endif

#Derive gfortran version from gcc version and check if >=4.6.0.

GVERSION := $(shell gcc -dumpversion)
GE40600 := $(shell expr `gcc -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40600)

ifeq ($(FC),gfortran)
  ifeq ($(GE40600), 0)
    $(error GCC/gfortran version $(GVERSION) detected. This program uses features of gfortran >= 4.6.0. Please update your toolchain)
  endif
else
  $(warning This code was optimized for gfortran >= 4.6.0 and has not been tested with other compilers)
endif

.SUFFIXES:
.SUFFIXES: .f90 .o .mod .out .test

PROGRAMS = maooam test_tl_ad

MODULES = util.mod params.mod inprod_analytic.mod aotensor_def.mod ic_def.mod integrator.mod stat.mod tensor.mod maooam_tl_ad.mod tl_ad_integrator.mod
MODULE_OBJECTS = $(MODULES:.mod=.o)

TEST_PROGRAMS = test_inprod_analytic test_aotensor
TESTS = $(TEST_PROGRAMS:=.test)

all: FCFLAGS = -llapack
all: $(PROGRAMS)

debug: FCFLAGS = -g -O0 -fbounds-check -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow
debug: all

test: FCFLAGS = -g -O0 -fbounds-check -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow
test: $(TESTS)

test_inprod_analytic.test: test_inprod_analytic
	./$< |sort  > tests/$<.out
	diff tests/$<.out tests/$<.ref && echo PASS: $@ || echo FAIL: $@

test_aotensor.test: test_aotensor 
	./$< |sort  > tests/$<.out
	diff tests/$<.out tests/$<.ref && echo PASS: $@ || echo FAIL: $@

%.mod: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.f90 $(MODULES)
	$(FC) $(FCFLAGS) -c $<

%: %.o $(MODULE_OBJECTS)
	$(FC) $(FCFLAGS) $^ -o $@

inprod_analytic.mod: params.mod
aotensor_def.mod: tensor.mod inprod_analytic.mod
integrator.mod: aotensor_def.mod
ic_def.mod: util.mod aotensor_def.mod
stat.mod: params.mod
maooam_tl_ad.mod: aotensor_def.mod
tl_ad_integrator.mod: maooam_tl_ad.mod

clean:
	$(RM) *.o *.mod $(PROGRAMS) $(TEST_PROGRAMS) tests/*.out

.PHONY: clean all test %.test
