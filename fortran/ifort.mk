# File containing the appropriate options for
# the ifort compiler

FC = ifort
FCFLAGS = -assume byterecl -O3
LDFLAGS = -mkl
DEBUG_FCFLAGS = -g -O0 -check all -traceback 

