# Makefile for transim2 input/output routines

include ../makefile.defs

SRCS = an_arith.c an_op_util.c variables.c an_dsp_op.c an_math_op.c \
	an_misc_op.c an_output.c an_optbl.c gnuplot.c report.c

CXXSRCS = init_clean.cc parser_functions.cc an_net_op.cc ../help.cc

OBJS = $(SRCS:.c=.o) $(CXXSRCS:.cc=.o) tr2_parse.tab.o lex.yy.o

TARGETLIB = ../lib/libio.a

$(TARGETLIB): $(TARGETLIB)($(OBJS)) 
	$(RANLIB) $(TARGETLIB)

tr2_parse.tab.o: tr2_parse.tab.c

tr2_parse.tab.c: tr2_parse.y
	$(YACC) -d tr2_parse.y

lex.yy.o: lex.yy.c

lex.yy.c: tr2_parse.lex tr2_parse.y
	$(LEX) -i tr2_parse.lex

.PHONY: clean
clean:
	-$(RM) $(OBJS) tr2_parse.tab.c tr2_parse.tab.h
	$(AR) d $(TARGETLIB) $(OBJS)

.PHONY: dep
dep:
	$(CC) $(CFLAGS) -MM $(SRCS) > .depend
	$(CXX) $(CXXFLAGS) -MM $(CXXSRCS) >> .depend

.depend: 
	echo > .depend

Makefile: ;

# Dependecies follow this line
-include .depend
