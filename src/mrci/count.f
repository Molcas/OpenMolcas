************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE COUNT(NINTGR,NSYM,NORB,MUL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NORB(*),MUL(8,8)
      DIMENSION NDPROD(8)
C COUNT TWO-ELECTRON INTEGRALS
C FIRST, COUNT NUMBER OF DENSITY PRODUCTS IN EACH SYMMETRY:
      NORBT=0
      DO 10 IS=1,NSYM
        NDPROD(IS)=0
        NORBT=NORBT+NORB(IS)
10    CONTINUE
      DO 30 IJS=1,NSYM
        ISUM=0
        DO 20 IS=1,NSYM
          JS=MUL(IS,IJS)
          IF(JS.GT.IS) GOTO 20
          ISUM=ISUM+NORB(IS)*NORB(JS)
20      CONTINUE
        NDPROD(IJS)=ISUM
30    CONTINUE
      NDPROD(1)=(NDPROD(1)+NORBT)/2
C THEN COUNT NUMBER OF TOTALLY SYMMETRIC PRODUCTS OF DENS-PRODUCTS:
      NINTGR=0
      DO 40 IJS=1,NSYM
        NINTGR=NINTGR+NDPROD(IJS)*(1+NDPROD(IJS))
40    CONTINUE
      NINTGR=NINTGR/2
      RETURN
      END
