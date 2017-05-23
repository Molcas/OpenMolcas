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
      SUBROUTINE IPO(IPOA,NVIR,MUL,NSYM,KLS,IFT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IPOA(*),NVIR(*),MUL(8,8)
      NSUM=0
      IF(IFT.LT.0) THEN
        DO 10 N=1,NSYM
          IPOA(N)=NSUM
          M=MUL(N,KLS)
          NSUM=NSUM+NVIR(N)*NVIR(M)
10      CONTINUE
      ELSE
        IF (KLS.EQ.1) THEN
          DO 20 N=1,NSYM
            IPOA(N)=NSUM
            NSUM=NSUM+(NVIR(N)*(NVIR(N)+1))/2
20        CONTINUE
        ELSE
          DO 30 N=1,NSYM
            IPOA(N)=NSUM
            M=MUL(N,KLS)
            IF(N.GT.M) THEN
              NSUM=NSUM+NVIR(N)*NVIR(M)
            END IF
30        CONTINUE
        END IF
      END IF
      IPOA(NSYM+1)=NSUM
      RETURN
      END
