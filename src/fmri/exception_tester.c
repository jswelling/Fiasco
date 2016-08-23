#include <stdio.h>
#include <fmri.h>
#include <fexceptions.h>

void callThings(int i)
{
  fprintf(stderr,"Hit callThings(%d)\n",i);
  if (i<6) callThings(i+1);
  else fex_raiseException(EXCEPTION_BASE,"Hello from level 6!");
  fprintf(stderr,"Returning from callThings(%d)\n",i);
}

void sub1()
{
  FEX_TRY ({
    callThings(0);
  }) 
    FEX_CATCH (EXCEPTION_IO, e,
    {
      fprintf(stderr,"I just got IO exception %d <%s>!\n",e->type,e->str);
      fprintf(stderr,"And now everything seems fine!\n");
      fex_raiseException(EXCEPTION_BASE,"Hello from sub1's exception hander!");
    });
  FEX_END_TRY;
}

int main(int argc, char* argv[])
{
  FEX_TRY
    ({
    fprintf(stderr,"I can do lots of stuff!\n");
    sub1();
    })
    FEX_CATCH(EXCEPTION_BASE, e,
    {
    fprintf(stderr,"I just got exception %s <%s>!\n",
	    fex_getExceptionTypeName(e), fex_getExceptionString(e));
    fprintf(stderr,"And now everything seems fine!\n");
    });
  FEX_END_TRY;


  fprintf(stderr,"All is now well!\n");
}
