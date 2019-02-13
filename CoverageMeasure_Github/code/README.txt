This code accompanies the publication

Sellentin & Starck (201x),"
	Debiasing inference with approximate covariance matrices and other unidentified biases", ArXiv: 1902.00709.

The code depends on the GSL (Gnu Scientific Library) and (in the future) OpenMP. The makefile assumes that both are installed in the compiler's/linker's standard search path, otherwise add flags -L and -I, as necessary.

Cover.h and Cover.cpp implement a base class. For standard uses, this base class is not intended to be modified.

The physical applications are specified in the *derived* classes, of which we include three examples: Example1.cpp, Example2.cpp, Example3.cpp. For your own application, either copy/paste/edit one of these classes, or add a new one.

Useful buzzwords to google include:
C++ inheritance, base class, derivation, c++ virtual functions, c++ std::vector, OpenMP manual.

Noteworthy wisdom: Depending on your data set, this code may push your university cluster's RAM and CPU usage to the limits. Never recompute anything which you already computed and please design your memory usage carefully.

Apart from this:

$ cd CoverageMeasure/code
$ make Cover.o                    #this compiles the base class
$ make Example1  && ./example1    #this compiles and runs Example1
$ make Example2  && ./example2    #this compiles and runs Example2
$ make Example3  && ./example3    #this compiles and runs Example3

There is no input parameter file. The codes need to be recompiled after changes.

