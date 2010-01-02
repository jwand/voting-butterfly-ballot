CC=g++

#
# add any your own library as necessary.
#
LIBS = 

#
# Two C source files you are asking to create for this project.
# You can add any additional files as necessary for this project.
#

SRCS.h = 
SRCS = fit_lqd1.cpp

#
# Object files will be created. If you add any additional source file
# , add here also.
#
OBJS = fit_lqd1.o


LDLIBS = -lm 

INCLUDEDIRS = -I/usr/lib/R/include -I/usr/local/lib/R/include
LDFLAGS = 

#
# header file used in this project.
#
CFLAGS = -O3 -ansi -Wunused  $(INCLUDEDIRS) $(DEBUGFLAG)

#
# Final program to be made. To run this project, simply type proj1.
#

PROGRAM =	fit_lqd.so

$(PROGRAM):	$(OBJS) $(RGENOUD.OBJS) $(LIBS)
		$(CC) -shared $(CFLAGS)  $(OBJS) $(RGENOUD.OBJS) $(LIBS) $(LDFLAGS) $(LDLIBS) -o $(PROGRAM)
		@echo type $(PROGRAM) to run.

fit_lqd1.o:	fit_lqd1.cpp $(SRCS.h)
		$(CC)  -c $(CFLAGS) fit_lqd1.cpp

# cleaning functions

clean:		
		rm -f $(OBJS) $(RGENOUD.OBJS) core

tclean:		
		rm -f $(OBJS) $(RGENOUD.OBJS) $(LIBS) $(PROGRAM) core 


