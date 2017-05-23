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
      SUBROUTINE CSCALE(INDX,INTSYM,C,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(*),INDX(*),INTSYM(*)

#include "SysDef.fh"

#include "mrci.fh"
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
CPAM96      JSYM(L)=UNPACK(INTSYM((L+9)/10),3*MOD(L-1,10)+1,3)+1
c      JSYM(L)=JSUNP(INTSYM,L)

c      DO 10 II1=IRC(3)+1,IRC(4)
       II1=IRC(3)+1
30       if(II1.gt.IRC(4)) goto 10
c        IF(JSYM(II1).NE.LSYM) GOTO 40
        IF(JSUNP(INTSYM,II1).NE.LSYM) GOTO 40
        NA=INDX(II1)
        MA=1
        if(NVIRT.lt.1) goto 620
c        DO 20 MA=1,NVIRT
720          C(NA+NDIAG(MA))=X*C(NA+NDIAG(MA))
c20      CONTINUE
        MA=MA+1
        if(MA.le.NVIRT) goto 720
620     continue
40      II1=II1+1
        GOTO 30
10    CONTINUE
      RETURN
      END
