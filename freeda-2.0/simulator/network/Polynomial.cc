#include "Polynomial.h"

Polynomial::Polynomial(const int& dimension, const int& no_of_coeff,double* coeff_values)
{
  contvoltvect_size = dimension ;
  coeffvect_size = no_of_coeff ;
  coefficient = coeff_values ;
}

inline int Polynomial::getDimension()
{
  return contvoltvect_size ;
}

AD Polynomial::findvalue(AD * contvoltvect)
{
  int i = 0 ;
  value = 0.0 ;

  if(coeffvect_size >= contvoltvect_size )
  {
    AD temp = 0 ;
    int int_count = 0 ;

    value = contvoltvect[0] ;   // first term

    for(i=0;i<contvoltvect_size;i++)
    {
      value += coefficient[i+1] * contvoltvect[i] ;       // second term
    }

    if(coeffvect_size == contvoltvect_size)
      return value ;

    k_count = contvoltvect_size + 1 ;
    value += subop(contvoltvect, 0) ;

    while(int_count != contvoltvect_size)
    {
      if (k_count > coeffvect_size)
      {
        value += value * contvoltvect[f_internal] ;
        return value ;
      }
      value += subop(contvoltvect,int_count) * contvoltvect[int_count];
      ++int_count ;
    }
  }
  return value ;
}

AD Polynomial::subop(AD * contvoltvect, int iterator)
{
  AD temp = 0.0 ;
  int count = 0 ;
  int internal = iterator ;

  f_internal = iterator ;

  while(count != contvoltvect_size)
  {
    if (k_count > coeffvect_size)
    {
      temp += (temp * contvoltvect[f_internal]) ;
      return temp ;
    }

    for(internal = iterator; internal < contvoltvect_size; internal++)
    {
      temp += coefficient[k_count] * contvoltvect[internal] ;
      ++k_count ;
      if (k_count > coeffvect_size)
      {
        temp += temp * contvoltvect[f_internal] ;
        return temp ;
      }
    }
    temp += temp * contvoltvect[f_internal] ;
    internal = f_internal + 1 ;
    f_internal = internal ;
    if(internal >= contvoltvect_size)
      break ;
    ++count ;
  }
  return temp ;
};
