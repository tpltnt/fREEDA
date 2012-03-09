# Elements makefile

PWD := $(shell pwd)
FRHOME=$(PWD)/../..

include ../makefile.defs
-include makefile.list

OBJS = $(SRCS:.cc=.o)

TARGETLIB = ../lib/libnetwork.a
YMATRIX_PATH = y/Ymatrix/src
NPORT_PATH = n/NPort/src
NPORTT_PATH = n/NPortT/src

$(TARGETLIB): $(TARGETLIB)($(OBJS)) element_headers.h
	$(RANLIB) $(TARGETLIB)

element_headers.h create_element.cc create_dummy_elem.cc: Makefile \
        mkheaders.sh makefile.list
	$(BASH) mkheaders.sh $(ELEMNAMES)

lex.Ymatrix_.c: $(YMATRIX_PATH)/Ymatrix.lex
	$(LEX) -B -PYmatrix_ -i $(YMATRIX_PATH)/Ymatrix.lex

lex.NPort_.c: $(NPORT_PATH)/NPort.lex
	$(LEX) -B -PNPort_ -i $(NPORT_PATH)/NPort.lex

lex.NPortT_.c: $(NPORTT_PATH)/NPortT.lex
	$(LEX) -B -PNPortT_ -i $(NPORTT_PATH)/NPortT.lex

.PHONY: clean
clean:
	-$(RM) $(OBJS)
	$(AR) d $(TARGETLIB) $(OBJS)
	-$(RM) element_headers.h create_element.cc create_dummy_elem.cc

.PHONY: dep
dep: lex.Ymatrix_.c lex.NPort_.c lex.NPortT_.c element_headers.h \
        create_element.cc create_dummy_elem.cc
	$(CXX) -MM $(CXXFLAGS) $(SRCS) > .depend

.depend: 
	echo >> .depend

makefile.list:
	$(CP) makefile.list.dist makefile.list

Makefile: ;


# Dependencies follow this line
-include .depend
