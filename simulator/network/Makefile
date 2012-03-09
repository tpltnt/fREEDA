# Network library makefile

include ../makefile.defs

SRCS = templates.cc \
	constants.cc \
	Circuit.cc \
	NetListItem.cc \
	CircuitManager.cc \
	Element.cc \
	Terminal.cc \
	GraphNode.cc \
	ElementManager.cc \
	ElementData.cc \
	ADInterface.cc \
	Polynomial.cc 


OBJS = $(SRCS:.cc=.o)

TARGETLIB = ../lib/libnetwork.a

$(TARGETLIB): .flagelem $(TARGETLIB)($(OBJS)) 
	$(RANLIB) $(TARGETLIB)

#
# We have the instances of all the templates here. So we compile
# without -fno-implicit-templates. In this way we avoid code
# duplication, at least by now. Besides, we do not generate tons of 
# debugging information for the template library.
#
templates.o: templates.cc
	$(CXX) $(CXXFLAGS) -c templates.cc

.PHONY: clean
clean:
	-$(RM) $(OBJS)
	$(AR) d $(TARGETLIB) $(OBJS)

.PHONY: dep
dep:
	$(CXX) -MM $(CXXFLAGS) $(SRCS) > .depend

.flagelem: ../elements/makefile.list
	$(MAKE) dep
	echo > .flagelem

.depend: 
	echo > .depend

Makefile: ;

# Dependencies follow this line
-include .depend

