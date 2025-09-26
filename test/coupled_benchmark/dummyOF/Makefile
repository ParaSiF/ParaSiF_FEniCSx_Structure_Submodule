#!/bin/bash

CPP	= mpic++
CC	= mpicc
CPPFLAGS	= -std=c++11 -O3
CFLAGS	= 
CPP_HEADER_PATH	= -I${PWD}/../../../../../coupling_lib/MUI
C_HEADER_PATH	=
CPP_LAB_PATH = 
C_LAB_PATH = 
CPP_LDLIBS = 
C_LDLIBS = -lstdc++ -lmpi_cxx -lmpi -lm
CPP_SCR = $(wildcard *.cpp)
C_SCR = $(wildcard *.c)
CPP_EXE = $(CPP_SCR:.cpp=)
C_EXE = $(C_SCR:.c=)


all: $(CPP_EXE) $(C_EXE)

% : %.cpp
	$(CPP) $(CPPFLAGS) $< -o $@ $(CPP_HEADER_PATH) $(CPP_LAB_PATH) $(CPP_LDLIBS)

% : %.c
	$(CC) $(CFLAGS) $< -o $@ $(C_HEADER_PATH) $(C_LAB_PATH) $(C_LDLIBS)


.PHONY: clean
	
clean:
	-rm -f $(CPP_EXE) $(C_EXE)