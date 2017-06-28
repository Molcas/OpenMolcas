#ifndef _HAVE_EXTRA_

      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n)
      REAL*8 d,a(np,np)
      END SUBROUTINE ludcmp

      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      END SUBROUTINE lubksb

#endif
