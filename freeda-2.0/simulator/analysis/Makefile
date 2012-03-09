# libanalysis.a makefile

include ../makefile.defs

SRCS =  AC.cc \
	Analysis.cc \
	anMain.cc \
	CircVector.cc \
	Euler.cc \
	FreqDomainSV.cc \
	FreqMNAM.cc \
	Nox_Interface.cc \
	slist.cc \
	SVTran.cc \
	SVTran2.cc \
	DC.cc \
	SVTr.cc \
	SVHB.cc \
	SVHB_phys.cc \
	TimeDomainSV.cc \
	TimeMNAM.cc \
	Tincidence.cc \
	Trapezoidal.cc 

CSRCS = iohb.c

OBJS = $(SRCS:.cc=.o) $(CSRCS:.c=.o)

TARGETLIB = ../lib/libanalysis.a

$(TARGETLIB): $(TARGETLIB)($(OBJS))
	$(RANLIB) $(TARGETLIB)

#
# We have the instances of all the templates here. So we compile
# without -fno-implicit-templates. In this way we avoid code
# duplication, at least by now. Besides, we do not generate tons of
# debugging information for the template library.
#
#templates.o: templates.cc
	#$(CXX) $(TEMPLATEFLAGS) -c templates.cc

.PHONY: clean
clean:
	-$(RM) $(OBJS)
	$(AR) d $(TARGETLIB) $(OBJS)

.PHONY: dep

dep:
	$(CXX) -MM $(CXXFLAGS) $(SRCS) > .depend

.depend:
	echo > .depend

Makefile: ;

# Dependencies follow this line
-include .depend

