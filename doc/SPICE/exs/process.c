#include <stdio.h>

main()
{
  char now;
  float sum, total, total_home, total_test, total_proj;
  int i;
  FILE *fin, *fout, *ffinal;

  fin = fopen("tmp","r");
  fout = fopen ("tmpout", "w");
   do
    {
       for (i=1; i<=3; i++)
	{
	  fscanf(fin, "%g", &sum); 
	  if (i==1) ;
	  else if (i==2) fprintf(fout, " %g ", sum);
	  else if (i==3) fprintf(fout, " %g ", sum);
	}
       fprintf(fout, "\n");
      }
    while (sum > -1);
  fclose(fin);
  fclose(fout);
}
	  
      
