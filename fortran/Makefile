
include ./gfortran.mk  # Rules for gfortran compiler
#include ./ifort.mk  # Rules for ifort compiler

SHELL = /bin/sh
DEBUG = false
#Adapt to host system
ifeq (,$(findstring Windows,$(OS)))
  RM = rm -f
else
  RM = del
endif



ifeq ($(FC),gfortran)
  #Derive gfortran version from gcc version and check if >=4.6.0.
  GVERSION := $(shell gcc -dumpversion)
  GE40600 := $(shell expr `gcc -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40600)
  ifeq ($(GE40600), 0)
    $(error GCC/gfortran version $(GVERSION) detected. This program uses features of gfortran >= 4.6.0. Please update your toolchain)
  endif
  CHG = sed -i --follow-symlink 's/USE IFPORT/\!USE IFPORT/g' *.f90; sed -i --follow-symlink 's/!*USE IFPORT/\!USE IFPORT/g' *.f90
else ifeq ($(FC),ifort)
  $(warning This code was optimized for ifort >= 14.0.2 and has not been tested with other compilers)
  CHG = sed -i --follow-symlink 's/\!*USE IFPORT/USE IFPORT/g' *.f90
else
  CHG = sed -i --follow-symlink 's/USE IFPORT/\!USE IFPORT/g' *.f90
  $(warning This code was optimized for gfortran >= 4.6.0 and and if >=14.0.2. It has not been tested with other compilers)
endif

.SUFFIXES:
.SUFFIXES: .f90 .o .mod .out .test

PROGRAMS = maooam_lyap maooam test_tl_ad maooam_lyap maooam_lyap_div

MODULES = util.mod params.mod inprod_analytic.mod aotensor_def.mod ic_def.mod rk2_integrator.mod stat.mod tensor.mod tl_ad_tensor.mod rk2_tl_ad_integrator.mod rk4_integrator.mod rk4_tl_ad_integrator.mod icdelta_def.mod lyap_vectors.mod rk4_integrator.mod lyap_stat.mod 
MODULE_OBJECTS = $(MODULES:.mod=.o)

TEST_PROGRAMS = test_inprod_analytic test_aotensor
TESTS = $(TEST_PROGRAMS:=.test)

all: pre-build $(PROGRAMS) 

pre-build: 
	$(CHG)

debug: FCFLAGS = $(DEBUG_FCFLAGS)
debug: all

test: FCFLAGS =$(DEBUG_FCFLAGS) 
test: $(TESTS)

test_inprod_analytic.test: test_inprod_analytic
	./$< |sort  > tests/$<.out
	diff tests/$<.out tests/$<.ref && echo PASS: $@ || echo FAIL: $@

test_aotensor.test: test_aotensor 
	./$< |sort  > tests/$<.out
	diff tests/$<.out tests/$<.ref && echo PASS: $@ || echo FAIL: $@

%.mod: %.f90
	$(FC) $(LDFLAGS) $(FCFLAGS) -c $<

%.o: %.f90 $(MODULES)
	$(FC) $(LDFLAGS) $(FCFLAGS) -c $<

%: %.o $(MODULE_OBJECTS)
	$(FC) $(FCFLAGS) $^ -o $@ $(LDFLAGS)

inprod_analytic.mod: params.mod
aotensor_def.mod: tensor.mod inprod_analytic.mod
rk2_integrator.mod: aotensor_def.mod
rk4_integrator.mod: aotensor_def.mod
ic_def.mod: util.mod aotensor_def.mod
stat.mod: params.mod
tl_ad_tensor.mod: aotensor_def.mod
rk2_tl_ad_integrator.mod: tl_ad_tensor.mod
lyap_vectors.mod: params.mod
rk4_tl_ad_integrator.mod: tl_ad_tensor.mod

clean:
	$(RM) *.o *.mod $(PROGRAMS) $(TEST_PROGRAMS) tests/*.out

.PHONY: clean all test %.test
