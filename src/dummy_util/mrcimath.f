#ifndef _HAVE_EXTRA_

      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n)
      REAL*8 d,a(np,np)
      CALL Unused_Real_Array(a)
      CALL Unused_Integer(n)
      CALL Unused_Integer(np)
      CALL Unused_Integer_Array(indx)
      CALL Unused_Real(d)
      END SUBROUTINE ludcmp

      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      CALL Unused_Real_Array(a)
      CALL Unused_Integer(n)
      CALL Unused_Integer(np)
      CALL Unused_Integer_Array(indx)
      CALL Unused_Real_Array(b)
      END SUBROUTINE lubksb

#endif
