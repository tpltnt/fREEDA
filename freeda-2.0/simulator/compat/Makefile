# Makefile for fREEDA input/output routines

include ../makefile.defs

SRCS = cx_dcomp.c ml.c st.c dsp.c equation.c in_interp.c

OBJS = $(SRCS:.c=.o) $(CXXSRCS:.cc=.o)

TARGETLIB = ../lib/libcompat.a

$(TARGETLIB): $(TARGETLIB)($(OBJS)) 
	$(RANLIB) $(TARGETLIB)

.PHONY: clean
clean:
	-$(RM) $(OBJS)
	$(AR) d $(TARGETLIB) $(OBJS)

.PHONY: dep
dep:
	$(CXX) $(CFLAGS) -MM $(SRCS) > .depend

.depend: 
	echo > .depend

Makefile: ;

# Dependecies follow this line
-include .depend
