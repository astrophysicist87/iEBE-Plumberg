SHELL=/bin/sh

SRCS= \
main.cpp \
Arsenal.cpp \
gauss_quadrature.cpp

HDRS= \
main.h \
Arsenal.h \
gauss_quadrature.h

MAKEFILE=makefile

COMMAND=interpolate_CF

OBJS= $(addsuffix .o, $(basename $(SRCS)))

CC= g++
CFLAGS=  -g
#CFLAGS=
WARNFLAGS= -ansi -pedantic -Wall -W \
   -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
   -Wwrite-strings -fshort-enums -fno-common 
LDFLAGS= -lgsl -lgslcblas 
LIBS= -L/sw/lib -I/sw/include

$(COMMAND): $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) $(CFLAGS) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJS)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
 
