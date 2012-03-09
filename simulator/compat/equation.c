/***********************************************************************
 * File Name:  equation.c
 *
 *
 *
 * To use as a subroutine cal expr()
 *	Expression evaluator
 *
 *
 * Major features:
 *   1. an unlimited number of variables
 *   2. supports basic trigonometric expressions
 *
 * Example: (> line entered by user, < system response)
 *      >  a = 3
 *	<  3
 *	>  a = a + a
 *	<  6
 *	>  3/2
 *	<  1.5
 *	> y = tan(2)
 *	< -2.18504
 *	>  (-2+3)*2
 *	< 2
 *	> --2
 *	< 2
 *
 * Parameters:
 *	string	  = text of expression terminated in null or '}'
 *	            text may begin with '{' but this is ignored.
 *	result	  = double precision result
 *	charsRead = count of number of characters read.
 *
 * Return codes:
 *	0	no error
 *	1	syntax error
 *	2	parameter not defined
 *	3	function syntax error
 *	4	generic function argument out of range
 *    101       ACOSH function argument out of range
 *    102       ATANH function argument out of range
 *    103       ACOS  function argument out of range
 *    104       ASIN  function argument out of range
 *    105       LOG   function argument out of range
 *    106       LOG10 function argument out of range
 *    107       SQRT  function argument out of range
 *    108       EXP   function argument out of range
 *    109       POWER (^ or **) function argument out of range
 *
 *
 *
 *
 * Note:
 *   The following constants can be used in expressions		  */
#define 	DTR	0.0174532925  /* degrees to radians conversion factor */
#define 	PI	3.141592654
/*
 *
 * Note:
 *	Standard expressions are evaluated and the following functions are
 *      supported.
 *
 *              OPERATOR		   SYNTAX */
typedef	enum	
{
  PLUS,			/* x+y			*/
  MINUS,			/* x-y			*/
  UNARY_PLUS,		/* +x			*/
  UNARY_MINUS,		/* -x			*/
  MULTIPLY,		/* x*y			*/
  DIVIDE,			/* y/x			*/	
  POW,			/* x^y  or   x**y	*/
  AND,			/* x&y			*/
  OR,			/* x|y			*/
  NOT,			/* !x			*/
  XOR,			/* ~x			*/
  /* Add Functions below this line and edit the following
   * Routines:
   *	evaluate()
   *	get_token()
   */
  SIN,			/* sin(x)	radians*/
  COS,			/* cos(x)	radians*/
  TAN,			/* tan(x)	radians*/
  ASIN,			/* asin(x)	radians*/
  ACOS,			/* acos(x)	radians*/
  ATAN,			/* atan(x)	radians*/
  SINH,			/* sinh(x)		*/
  COSH,			/* cosh(x)		*/
  TANH,			/* tanh(x)		*/
  EXP,			/* exp(x)	i.e. e^x*/
  ASINH,			/* asinh(x)		*/
  ACOSH,			/* acosh(x)		*/
  ATANH,			/* atanh(x)		*/
  ABS,			/* abs(x)		*/
  SQRT,			/* sqrt(x)		*/
  LOG,			/* log(x)		*/
  LOG10,			/* log10(x)		*/
  SIGN,			/* sign(x)		*/
  DUPLICATE,		/* dup(x)  duplicates x on stack*/
  /* Add Functions above this line. */
  LEFT_PAREN,
  RIGHT_PAREN,
}	Operator;


#include	<stdio.h>
#include	<stdlib.h>
#include	<ctype.h>
#include	<math.h>
#include	<float.h>
#include        <string.h>
#include	"standard.h"
#include	"mltypes.h"
#include	"generic.h"
#include	"st.h"
#include	"equation.h"

#define		LN_MAXREAL	87.5
#define		MAXREAL		FLT_MAX
#define		MINREAL		FLT_MIN


typedef	double		Operand;

static	Operand value, value1;

typedef	enum	
{
  OPERATOR,
  FUNCTION,
  OPERAND,
  DELIMITER,
  END_OF_EXPR
}	TType;

typedef	enum	
{
  ASSIGN
}	Delimiter;

typedef	struct	
{
  TType	type;
  union	
  {
    Operator operator;
    Operand operand;
    Delimiter delimiter;
  } val;
}	Token;


static void	push_operator(), push_operand(), get_text(), evaluate(),
  skip_space(), fatal_error();
static char	get_byte(), peek_byte();
static int	precedence();
static Operator	pop_operator();
static Operand	pop_operand();
static Token	get_token();
static int	isdigit_or_dot();

#define		TOStype		operatorType_stack[operator_sp]
#define		TOSisFUNCTION	(operator_stack[operator_sp] < 0) ? 1 : 0
#define		TOS		operator_stack[operator_sp]
#define		EMPTY		-1

#define		push_function(x)	push_operator(-x)

static char *input_buf;

static Operand	operand_stack[STACK_SIZE];
static int	operand_sp;

static Operator	operator_stack[STACK_SIZE];
static int	operator_sp;

static int	err_flag;
static int	charCount, oldCharCount;



/*
 * strdiff
 * caseless strcmp function
 */
int strdiff(a,b) /* returns 0 if string a == string b */
     char *a, *b;
{
  while ((*a++ & 0X5F) == (*b++ & 0X5F)) if (!(*a)) return(*b);
  return(1);
}

int expr(string, result, charsRead)
     char *string;  /* null terminated INPUT */
     double *result;  /* OUTPUT */
     int *charsRead;  /* OUTPUT */
{

  Token	token;
  /* Operand	result; */
  TType	history;

  //printf("+++ expr() enter\n");
  //printf("+++ expr() string = %s\n",string);
#ifdef DEBUG
  Eq_DumpSymbolTable(stdout, capOptionsT_P);
#endif

  charCount = oldCharCount = 0;
  operator_sp = operand_sp = EMPTY;
  history = OPERATOR;
  err_flag = 0;

  input_buf = string;

  do	
  {
    token = get_token();
    //printf("+++ expr() token.type =%d\n",token.type);
    if(err_flag) goto ERROR;
    switch (token.type)
    {
      /*        case DELIMITER:*/
      /*switch (token.val.delimiter)*/
      /*{*/
      /*case LEFT_BRACE:*/
      /*if( (operator_sp != EMPTY)*/
      /*&& (operand_sp != EMPTY))*/
      /*{*/
      /*err_flag = 1;*/
      /**result = 0.0;*/
      /**charsRead = charCount;*/
      /*return(err_flag);*/
      /*}*/
      /*break;  */
      /*case RIGHT_BRACE:*/
      /*token.type = END_OF_EXPR;*/
      /*break;*/
      /*}*/
      /*break;*/
      case OPERATOR:
        if (history == OPERATOR)
        {
          switch (token.val.operator)
          {
            case MINUS:
              token.val.operator =
                UNARY_MINUS;
              break;
            case PLUS:
              token.val.operator =
                UNARY_PLUS;
              break;
            default:
              ;
          }
        }
        switch (token.val.operator)
        {
          case LEFT_PAREN:
          case UNARY_PLUS:
          case UNARY_MINUS:
          case NOT:
            push_operator(token.val.operator);
            if(err_flag) goto ERROR;
            break;
          case RIGHT_PAREN:
            while (TOS != LEFT_PAREN)
            {
              evaluate();
              if(err_flag) goto ERROR;
            }
            pop_operator();
            if(err_flag) goto ERROR;
            token.type = OPERAND;
            if(operator_sp != EMPTY)
              if(TOSisFUNCTION)
              {
                evaluate();
                if(err_flag) goto ERROR;
              }
            break;
          default:
            while ((operator_sp != EMPTY) &&
                (precedence(token.val.operator)
                 <= precedence(TOS)))
            {
              evaluate();
              if(err_flag) goto ERROR;
            }
            push_operator(token.val.operator);
            if(err_flag) goto ERROR;
        }
        break;
      case FUNCTION:
        //printf("+++expr Function\n");
        push_function(token.val.operator);
        if(err_flag) goto ERROR;
        break;
      case OPERAND:
        push_operand(token.val.operand);
        if(err_flag) goto ERROR;
        break;
      default:
        ;
    }
    history = token.type;
    if (err_flag)
    {
      /* REPORT(WARNING_FATAL_INTERNAL,"Error in expression");
      */
      fatal_error("Error in expression.");
      return(0);
    }
  } while (token.type != END_OF_EXPR);

  while (operator_sp != EMPTY)
  {
    evaluate();
    if(err_flag) goto ERROR;
  }
  *result = pop_operand();
  if (operand_sp != EMPTY) err_flag = 1;
  if(err_flag) goto ERROR;
  *charsRead = charCount;
#ifdef DEBUG
  printf("equation: result = %12.2g\n", *result);
#endif
  return(err_flag);

ERROR:
  {
    char err_string[200];
    int length;
    *result = 0.0;		
    *charsRead = oldCharCount;	
    switch(err_flag)
    {
      case 1:
        sprintf(err_string,
            "Expression evaluation error: syntax error in\n   %s",string);
        break;
      case 2:
        sprintf(err_string,
            "Expression evaluation error: undefined parameter in\n   %s",string);
        break;
      case 3:
        sprintf(err_string,
            "Expression evaluation error: function syntax error in\n   %s" ,string);
        break;
      case 4:
        sprintf(err_string,
            "\n  Expression evaluation error: function argument out of range in\n   %s"
            ,string);
        break;
      case 101:
      case 102:
      case 103:
      case 104:
      case 105:
      case 106:
      case 107:
      case 108:
      case 109:
        sprintf(err_string, "\n  Expression evaluation error: ");
        length = strlen(err_string);
        switch(err_flag)
        {
          case 101:
            sprintf(&err_string[length],"ACOSH(%g)",value);
            break;
          case 102:
            sprintf(&err_string[length],"ATANH(%g)",value);
            break;
          case 103:
            sprintf(&err_string[length],"ACOS(%g)",value);
            break;
          case 104:
            sprintf(&err_string[length],"ASIN(%g)",value);
            break;
          case 105:
            sprintf(&err_string[length],"LOG(%g)",value);
            break;
          case 106:
            sprintf(&err_string[length],"LOG10(%g)",value);
            break;
          case 107:
            sprintf(&err_string[length],"SQRT(%g)",value);
            break;
          case 108:
            sprintf(&err_string[length],"EXP(%g)",value);
            break;
          case 109:
            sprintf(&err_string[length],
                "\n   POWER (%g^%g or %g**%g)",value1,value,value1,value);
            break;
        }
        length = strlen(err_string);
        sprintf(&err_string[length],
            " function argument out of range in\n   %s\n",string);
    }
    fatal_error(err_string);	
  }
  //printf("+++ expr() leave\n");
  return(err_flag);
}

static void push_operator(x)
     Operator x;
{
  if (operator_sp >= STACK_SIZE-1)
    err_flag = 1;
  else
    operator_stack[++operator_sp] = x;
}


static void push_operand(x)
     Operand x;
{
  if (operand_sp >= STACK_SIZE-1)
    err_flag = 1;
  else
    operand_stack[++operand_sp] = x;
}

static Operator pop_operator()
{
  if (operator_sp < 0)
  {
      err_flag = 1;
      return(0);
  }
  else
    return(abs(operator_stack[operator_sp--]));
}

static Operand pop_operand()
{
  if (operand_sp < 0)
  {
      err_flag = 1;
      return(0);
  }
  else
    return(operand_stack[operand_sp--]);
}

static int precedence(operator)
     Operator operator;
{
  switch (operator)
  {
    case LEFT_PAREN:
    case RIGHT_PAREN:
      return(0);
    case XOR:
      return(2);
    case OR:
      return(3);
    case AND:
      return(4);
    case PLUS:
    case MINUS:
      return(5);
    case MULTIPLY:
    case DIVIDE:
      return(6);
    case POW:
      return(7);
    case UNARY_PLUS:
    case UNARY_MINUS:
    case NOT:
      return(8);
    default: /* Function assumed */
      return(1);
  }
}

static void evaluate()
{
  int  ivalue, ivalue1;
  switch(pop_operator())
  {
    case XOR:
      ivalue = (pop_operand() == 0.0) ? 0 : 1;
      ivalue1 = (pop_operand() == 0.0) ? 0 : 1;
      value = (ivalue && !ivalue1) || (!ivalue && ivalue1);
      push_operand(value);
      break;
    case OR:
      ivalue = (pop_operand() == 0.0) ? 0 : 1;
      ivalue1 = (pop_operand() == 0.0) ? 0 : 1;
      value = (ivalue || ivalue1);
      push_operand(value);
      break;
    case AND:
      ivalue = (pop_operand() == 0.0) ? 0 : 1;
      ivalue1 = (pop_operand() == 0.0) ? 0 : 1;
      value = (ivalue && ivalue1);
      push_operand(value);
      break;
    case PLUS:
      push_operand(pop_operand() + pop_operand());
      break;
    case MINUS:
      value = pop_operand();
      push_operand(pop_operand() - value);
      break;
    case MULTIPLY:
      push_operand(pop_operand() * pop_operand());
      break;
    case DIVIDE:
      value = pop_operand();
      push_operand(pop_operand() / value);
      break;
    case UNARY_PLUS:
      push_operand(pop_operand());
      break;
    case UNARY_MINUS:
      push_operand(-pop_operand());
      break;
    case POW:
      value = pop_operand();
      value1 = pop_operand();
      if(value1<0)
      {
        ivalue = floor(value);
        if(fabs(ivalue-value) < MINREAL)
        {
          value = pow(-value1, value);
          if(ivalue%2) value = -value;
          push_operand(value);
        }
        else
          err_flag=109;
      }
      else
        push_operand(pow(value1, value));
      break;
      /* Functions are below this line. */
    case SIN:
      push_operand(sin(pop_operand()));
      break;
    case COS:
      push_operand(cos(pop_operand()));
      break;
    case TAN:
      push_operand(tan(pop_operand()));
      break;
    case ASIN:
      value = pop_operand();
      if(value > 1.0)
        err_flag=104;
      else
        push_operand(asin(value));
      break;
    case ACOS:
      value = pop_operand();
      if(value > 1.0)
        err_flag=103;
      else
        push_operand(acos(value));
      break;
    case ATAN:
      push_operand(atan(pop_operand()));
      break;
    case SINH:
      push_operand(sinh(pop_operand()));
      break;
    case COSH:
      push_operand(cosh(pop_operand()));
      break;
    case TANH:
      push_operand(tanh(pop_operand()));
      break;
    case ASINH:
      value = pop_operand();
      push_operand(
		   log(value + sqrt(value * value + 1)));
      break;
    case ACOSH:
      value = pop_operand();
      if(value <= 1.0)
        err_flag = 101;
      else
        push_operand(log(value + sqrt(value * value - 1)));
      break;
    case ATANH:
      value = pop_operand();
      if(fabs(value) >= 1.0)
        err_flag = 102;
      else
        push_operand(0.5 * log((1.0 + value) / (1.0-value)));
      break;
    case ABS:
      push_operand(fabs(pop_operand()));
      break;
    case EXP:
      value = pop_operand();
      if(value > LN_MAXREAL)
        err_flag=108;
      else
        push_operand(exp(value));
      break;
    case SQRT:
      value = pop_operand();
      if(value < 0)
        err_flag = 107;
      else
        push_operand(sqrt(value));
      break;
    case LOG:
      value = pop_operand();
      if(value > MAXREAL || value < MINREAL)
        err_flag = 105;
      else
        push_operand(log(value));
      break;
    case LOG10:
      value = pop_operand();
      if(value > MAXREAL || value < MINREAL)
        err_flag = 106;
      else
        push_operand(log10(value));
      break;
    case SIGN:
      push_operand(pop_operand() >= 0.0 ? 1.0 : -1.0);
      break;
    case DUPLICATE:
      value = pop_operand();
      push_operand(value);
      push_operand(value);
      break;
      /* Functions are above this line. */
    case NOT:
      value = (pop_operand() == 0.0) ? 1 : 0;
      push_operand(value);
      break;
    case LEFT_PAREN:
    case RIGHT_PAREN:
      err_flag = 1;
  }
}

static Token get_token()
{
  Token token;
  char text[MAX_STRING_LEN];
  double atof();

  //printf("+++ get_token enter\n");
  //printf("+++ get_token charCount = %d\n",charCount);
  oldCharCount = charCount;
  get_text(text);
  //printf("+++ get_token text = %s\n",text);
  switch (text[0])
  {
    case '\0':
      token.type = END_OF_EXPR;
      break;
    case '^':
      token.type = OPERATOR;
      token.val.operator = POW;
      break;
    case '~':
      token.type = OPERATOR;
      token.val.operator = XOR;
      break;
    case '!':
      token.type = OPERATOR;
      token.val.operator = NOT;
      break;
    case '|':
      token.type = OPERATOR;
      token.val.operator = OR;
      break;
    case '&':
      token.type = OPERATOR;
      token.val.operator = AND;
      break;
    case '+':
      token.type = OPERATOR;
      token.val.operator = PLUS;
      break;
    case '-':
      token.type = OPERATOR;
      token.val.operator = MINUS;
      break;
    case '*':
      token.type = OPERATOR;
      token.val.operator = MULTIPLY;
      break;
    case '/':
      token.type = OPERATOR;
      token.val.operator = DIVIDE;
      break;
    case '(':
      token.type = OPERATOR;
      token.val.operator = LEFT_PAREN;
      break;
    case ')':
      token.type = OPERATOR;
      token.val.operator = RIGHT_PAREN;
      break;
    case '{':
      token.type = OPERATOR;
      token.val.operator = LEFT_PAREN;
      break;
    case '}':
      token.type = OPERATOR;
      token.val.operator = RIGHT_PAREN;
      break;
    case '=':
      token.type = DELIMITER;
      token.val.delimiter = ASSIGN;
      break;
    default:
      if (isdigit_or_dot(text[0]))
      {
        token.type = OPERAND;
        token.val.operand = atof(text);
      }
      /* Functions are below this line. */
      else
        if(!strdiff(text, "PI"))  {
          token.type = OPERAND;
          token.val.operand = PI;
        }
        else if(!strdiff(text, "DTR"))  {
          token.type = OPERAND;
          token.val.operand = DTR;
        }
        else if(!strdiff(text, "SIN"))  {
          token.type = FUNCTION;
          token.val.operator = SIN;
        }
        else if(!strdiff(text, "COS"))  {
          token.type = FUNCTION;
          token.val.operator = COS;
        }
        else if(!strdiff(text, "TAN"))  {
          token.type = FUNCTION;
          token.val.operator = TAN;
        }
        else if(!strdiff(text, "ASIN"))  {
          token.type = FUNCTION;
          token.val.operator = ASIN;
        }
        else if(!strdiff(text, "ACOS"))  {
          token.type = FUNCTION;
          token.val.operator = ACOS;
        }
        else if(!strdiff(text, "ATAN"))  {
          token.type = FUNCTION;
          token.val.operator = ATAN;
        }
        else if(!strdiff(text, "SINH"))  {
          token.type = FUNCTION;
          token.val.operator = SINH;
        }
        else if(!strdiff(text, "COSH"))  {
          token.type = FUNCTION;
          token.val.operator = COSH;
        }
        else if(!strdiff(text, "TANH"))  {
          token.type = FUNCTION;
          token.val.operator = TANH;
        }
        else if(!strdiff(text, "ASINH"))  {
          token.type = FUNCTION;
          token.val.operator = ASINH;
        }
        else if(!strdiff(text, "ACOSH"))  {
          token.type = FUNCTION;
          token.val.operator = ACOSH;
        }
        else if(!strdiff(text, "ATANH"))  {
          token.type = FUNCTION;
          token.val.operator = ATANH;
        }
        else if(!strdiff(text, "ABS"))  {
          token.type = FUNCTION;
          token.val.operator = ABS;
        }
        else if(!strdiff(text, "SQRT"))  {
          token.type = FUNCTION;
          token.val.operator = SQRT;
        }
        else if(!strdiff(text, "EXP"))  {
          token.type = FUNCTION;
          token.val.operator = EXP;
        }
        else if(!strdiff(text, "LOG"))  {
          token.type = FUNCTION;
          token.val.operator = LOG;
        }
        else if(!strdiff(text, "LOG10"))  {
          token.type = FUNCTION;
          token.val.operator = LOG10;
        }
        else if(!strdiff(text, "SIGN"))  {
          token.type = FUNCTION;
          token.val.operator = SIGN;
        }
        else if(!strdiff(text, "DUP"))  {
          token.type = FUNCTION;
          token.val.operator = DUPLICATE;
        }
      /* Functions are above this line. */
        else { /* This should be a variable */
          if(St_GetSymAsDouble(capOptionsT_P, text,
                &token.val.operand) != 
              ST_SYM_FOUND)
          {
            err_flag = 2;
            return(token);
          }
          token.type = OPERAND;
#ifdef DEBUG
          printf("equation: variable usage, variable %s = %g\n",
              text, token.val.operand);
#endif
        }
  }
  if(token.type == FUNCTION)
    if(peek_byte() != '(')
      err_flag = 3;
  return(token);
}

static int ispunctuation(byte)
     char byte;
{
  if(   byte == '^'
	|| byte == '~'
	|| byte == '='
	|| byte == '!'
	|| byte == '|'
	|| byte == '&'
	|| byte == '+'
	|| byte == '-'
	|| byte == '*'
	|| byte == '/'
	|| byte == ','
	|| byte == ' '
	|| byte == '\t'
	|| byte == '\n'
	|| byte == '\0'
	|| byte == '('
	|| byte == ')') return(1);
  else return(0);
}

static void get_text(text)
     char *text;
{
  /* char get_byte(), peek_byte(); */
  skip_space();
  if ((*text = get_byte()) != '\0')
  {
    //printf("+++get_text() A: *text = %c\n",*text);
    if (*text == '*' && peek_byte() == '*')
    {
      *text++ = '^';
      get_byte();
    }
    else if (isdigit_or_dot(*text))
    {
      text++;
      while (isdigit_or_dot(peek_byte()))
        *text++ = get_byte();
      if (peek_byte() == 'e' || peek_byte() == 'E' ||
          peek_byte() == 'd' || peek_byte() == 'D')
      {
        get_byte();
        *text++ = 'e';
        *text++ = get_byte();	/* may be + or - */
        while (isdigit(peek_byte()))
          *text++ = get_byte();
      }
    }
    else if (isalpha(*text))
    {
      text++;
      while (!ispunctuation(peek_byte()))
      {
        if((*text = get_byte()) == '\0') return;
        // printf("+++get_text() B: *text = %c\n",*text);
        text++;
      }
    }
    else text++;
    *text = '\0';
  }
  //printf("+++get_text(): leave\n");
}

static int isdigit_or_dot(c)
     char c;
{
  if (c == '.' || isdigit(c)) return(1);
  else return(0);
}


static void skip_space()
{
  while (1)
  {
    switch (*input_buf)
    {
      case ' ':
      case '\t':
      case '\n':
        ++input_buf;
        break;
      default:
        return;
    }
  }
}

static char get_byte()
{
  if (*input_buf == '\0') return('\0');
  else
  {
    charCount++;
    return(*input_buf++);
  }
}

static char peek_byte()
{
  return(*input_buf);
}

void init_equation()
{

}

void free_equation()
{
  St_DelTable(capOptionsT_P);
}

int equation(string, result)
  char *string;
  double *result;
{
  int position;
  int status;
  int assign, assignCharCount;
  char assignName[MAX_STRING_LEN];

  /* Check to see if a string followed by assignment */
  {
    char text[MAX_STRING_LEN];
    input_buf = string;
    assignCharCount = charCount = oldCharCount = 0;
    get_text(text);
    if(isalpha(text[0]))
    {
      strncpy(assignName, text, MAX_STRING_LEN-1);
      get_text(text);
      if(text[0] == '=')
      {
        assign = 1;
        string = input_buf;
        assignCharCount = charCount;
      }
      else
        assign = 0;
    }
    else
      assign = 0;
  }

  if( (status=expr(string, result, &position)) )
  {
    char error_string[50];
    position += assignCharCount;
#ifdef DEBUG
    printf("charsRead = %d\n", position);
#endif
    switch(status)
    {
      case 1:
        sprintf(error_string,"Syntax error in expression.");
        break;
      case 2:
        sprintf(error_string,"Parameter not defined in expression.");;
        break;
      case 3:
        sprintf(error_string,"Function syntax error in expression.");;
        break;
      default:
        sprintf(error_string,"Error in expression.");
    }
    fatal_error(error_string);
    *result = 0;
  }

  if(assign)
  {
    generic_t gval;
    gval.d = *result;
    St_ReplSym(capOptionsT_P, assignName, GEN_TYPE_DOUBLE, gval);
  }
  return(status);
}


static void fatal_error(text)
  char *text;
{
  /*	PrintErr(text);
  */
  printf("%s\n",text);
}


/*
 * Dump a symbol table
 * (for debugging or utility purposes)
 */
//#ifdef DEBUG
void Eq_DumpSymbolTable(FILE *output_F, st_Table_Pt st_P)
{
  st_Entry_Pt	entry_P;
  int		bucket, i, size;
  usageId_Pt	usage_P;
  int		index;

  if(!st_P)
  {
    fprintf(output_F, "  Table ptr is NULL.\n");
    return;
  }

  fprintf(output_F, "  '%s' table, %d parameters defined\n", st_P->name,
      st_P->entryCt);

  for (bucket = 0; bucket < st_P->size; bucket++)
  {
    entry_P = st_P->tab_A[bucket];
    while (entry_P)
    {
      fprintf(output_F, "  %d: \t", bucket);

      switch(entry_P->type)
      {
        case GEN_TYPE_INT:
          fprintf(output_F, "  '%s' \t= %d \t(int)\n", entry_P->name,
              entry_P->gval.i);
          break;
        case GEN_TYPE_LONG:
          fprintf(output_F, "  '%s' \t= %ld \t(long)\n", entry_P->name,
              entry_P->gval.l);
          break;
        case GEN_TYPE_FLOAT:
          fprintf(output_F, "  '%s' \t= %g \t(float)\n", entry_P->name,
              entry_P->gval.f);
          break;
        case GEN_TYPE_DOUBLE:
          fprintf(output_F, "  '%s' \t= %g \t(double)\n", entry_P->name,
              entry_P->gval.d);
          break;
        case GEN_TYPE_DOUBLE_VECTOR:
          size = entry_P->gval.dv.size;
          fprintf(output_F, "  '%s' \t \t(double vector) size = %d \n",
              entry_P->name, size);
          for(i = 0; i < size; i++)
            fprintf(output_F, "\t %g\n", entry_P->gval.dv.v[i]);
          break;
        case GEN_TYPE_CHAR:
          fprintf(output_F, "  '%s' \t= '%c'/%d \t(char)\n", entry_P->name,
              entry_P->gval.c, entry_P->gval.c);
          break;
        case GEN_TYPE_STRING:
          if(entry_P->gval.s)
            fprintf(output_F, "  '%s' \t= '%s' \t(string)\n", 
                entry_P->name, entry_P->gval.s);
          else
            fprintf(output_F, "  '%s' \t= NULL \t(string)\n", 
                entry_P->name);
          break;
        case GEN_TYPE_STRING_A:
          if(!entry_P->gval.s_A)
            fprintf(output_F, "  '%s' \t= NULL \t(string array)\n",
                entry_P->name);
          else
          {
            char	**s_A;
            fprintf(output_F, "  '%s' is string array:\n",
                entry_P->name);
            for (s_A = entry_P->gval.s_A; *s_A; s_A++)
              fprintf(output_F, "\t\t '%s'\n", *s_A);
          }
          break;
        case GEN_TYPE_VOID_PTR:
          //	      if(entry_P->gval.v_P)
          //		fprintf(output_F, "  '%s' \tis user type \t(void ptr)\n",
          //			entry_P->name);
          //    break;
          //      }
          fprintf(output_F, "  '%s' \tis NULL user type \t(void ptr)\n",
              entry_P->name);
          break;
      default:
          fprintf(output_F, "  '%s' \tis unknown type\n", entry_P->name);
          break;
    }

    /* Get next entry. */
    entry_P = entry_P->next_P;
  }
}
return;
}
//#endif

