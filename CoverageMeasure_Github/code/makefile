# Makefile for coverage calibration of posteriors,
# or other methods which produce credibility contours.
#
# Ln(a) Sellentin
# Sterrewacht, Universiteit Leiden
# 2018



#LINKER FLAGS AND FLAGS TO BE TOUCHED~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# GSL_LDPATH: The path to where the files libgsl.so and/or libgsl.a lie. 
GSL_LDPATH = #-L/usr/local/lib -L/ ??? /gsl-1.16/.libs/

# GSL_INCDIR: The path to the folder that contains the folder
# `gsl' that contains e.g. the file gsl_math.h
GSL_INCDIR = #-I/usr/local/include/ -I/ ??? /gsl-1.16/include/

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Settings for C++ compiler and linker
CC = g++
CCFLAGS = -std=c++98 -O2  #$(GSL_INCDIR) -pg -no-pie
#LINKER FLAGS
LDFLAGS	 = $(GSL_LDPATH) -lgsl -lgslcblas 

U = utils/


################### OBJECTS ####################################################

Example1: Example1.cpp Cover.o Matrices.o
	$(CC) $(CCFLAGS) -c -o Example1.o Example1.cpp  \
	&& $(CC) $(CCFLAGS) Example1.o Cover.o $(U)Matrices.o $(LDFLAGS) -o example1

MeasureVariance: MeasureVariance.cpp Cover.o Matrices.o
	$(CC) $(CCFLAGS) -c -o MeasureVariance.o MeasureVariance.cpp  \
	&& $(CC) $(CCFLAGS) MeasureVariance.o Cover.o $(U)Matrices.o $(LDFLAGS) -o measurevariance
	 
Example2: Example2.cpp Cover.o Matrices.o
	$(CC) $(CCFLAGS) -c -o Example2.o Example2.cpp  \
	&& $(CC) $(CCFLAGS) Example2.o Cover.o $(U)Matrices.o $(LDFLAGS) -o example2

Example3: Example3.cpp Cover.o Matrices.o
	$(CC) $(CCFLAGS) -c -o Example3.o Example3.cpp  \
	&& $(CC) $(CCFLAGS) Example3.o Cover.o $(U)Matrices.o $(LDFLAGS) -o example3

Cover.o: Cover.cpp Cover.h Matrices.o
	$(CC) $(CCFLAGS) -c -o Cover.o Cover.cpp

Matrices.o: $(U)Matrices.cpp $(U)Matrices.h
	$(CC) $(CCFLAGS) -c -o $(U)Matrices.o $(U)Matrices.cpp


clean:
	rm *.o && rm utils/*.o


