# Freeda makefile

-include makefile.defs

SUBDIRS = elements network analysis compat inout

SRCS = freeda.cc help.cc 

OBJS = $(SRCS:.cc=.o)

# This is intended to recompile freeda.cc and 
# update the version every time.
.INTERMEDIATE: freeda.o

freeda: libs $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) -o freeda $(LDFLAGS)

all: dep freeda

.PHONY: libs
libs:
	if [ ! -d lib ]; then \
	  $(MKDIR) lib;\
	fi
	for DIRECT in $(SUBDIRS); do \
	  cd $$DIRECT && $(MAKE) && cd .. ; \
	done

makefile.defs: #config.h
	if [ ! -e makefile.defs ]; then \
	  $(CP) makefile.defs.sample makefile.defs; \
	fi
	@echo
	@echo ----------------------------------------------------------
	@echo
	@echo Run "make dep; make" to compile for the first time or when
	@echo adding elements.
	@echo
	@echo ----------------------------------------------------------
	@echo
	@false

#config.h:
	#$(CP) config.h.sample config.h

.PHONY: clean
clean:
	-$(RM) $(OBJS) freeda
	for DIRECT in $(SUBDIRS); do \
	  cd $$DIRECT && $(MAKE) clean && cd .. ; \
	done

.PHONY: cleanall
cleanall: clean
	-$(RM) lib/*.a

.PHONY: dep
dep: 
	$(CXX) $(OPTFLAGS) -MM $(SRCS) > .depend
	for DIRECT in $(SUBDIRS); do \
	  cd $$DIRECT && $(MAKE) dep && cd .. ; \
	done

.depend:
	echo > .depend

Makefile: ;


# Dependecies follow this line
-include .depend
