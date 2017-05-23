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
      SUBROUTINE GETINCN_RASSCFS(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,INTLST,NSMOB,ICOUL)
*
* Obtain integrals
*
*     ICOUL = 0 :
*                  XINT(IK,JL) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 1 :
*                  XINT(IJ,KL) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 2 :  XINT(IL,JK) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
* Storing for ICOUL = 1 not working if IKSM or JLSM .ne. 0
*
*
* Version for integrals stored in INTLST
*
* If type equals zero, all integrals of given type are fetched
* ( added aug8, 98)
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.fh"
#include "orbinp.fh"
*. Integral list
      Real * 8 Intlst(*)
*.Output
      DIMENSION XINT(*)
*. Local scratch
C-jwk-cleanup      DIMENSION IJARR(MXPORB)
*
      IF(ITP.GE.1) THEN
        iOrb=NOBPTS(ITP,ISM)
      ELSE
        IORB = NTOOBS(ISM)
      END IF
*
      IF(JTP.GE.1) THEN
        jOrb=NOBPTS(JTP,JSM)
      ELSE
        JORB = NTOOBS(JSM)
      END IF
*
      IF(KTP.GE.1) THEN
        kOrb=NOBPTS(KTP,KSM)
      ELSE
        KORB = NTOOBS(KSM)
      END IF
*
      IF(LTP.GE.1) THEN
        lOrb=NOBPTS(LTP,LSM)
      ELSE
        LORB = NTOOBS(LSM)
      END IF
*. Offsets relative to start of all orbitals, symmetry ordered
      IOFF = IBSO(ISM)
      DO IITP = 1, ITP -1
        IOFF = IOFF + NOBPTS(IITP,ISM)
      END DO
*
      JOFF = IBSO(JSM)
      DO JJTP = 1, JTP -1
        JOFF = JOFF + NOBPTS(JJTP,JSM)
      END DO
*
      KOFF = IBSO(KSM)
      DO KKTP = 1, KTP -1
        KOFF = KOFF + NOBPTS(KKTP,KSM)
      END DO
*
      LOFF = IBSO(LSM)
      DO LLTP = 1, LTP -1
        LOFF = LOFF + NOBPTS(LLTP,LSM)
      END DO

*
*     Collect Coulomb terms
*
*
      iInt=0
      Do l=lOff,lOff+lOrb-1
        jMin=jOff
        If ( JLSM.ne.0 ) jMin=l
        Do j=jMin,jOff+jOrb-1
          Do k=kOff,kOff+kOrb-1
            iMin = iOff
            If(IKSM.ne.0) iMin = k
            IF(ICOUL.EQ.1)  THEN
*. Address before integral (1,j!k,l)
                IINT = (L-LOFF)*Jorb*Korb*Iorb
     &               + (K-KOFF)*Jorb*Iorb
     &               + (J-JOFF)*Iorb
            ELSE IF (ICOUL.EQ.2) THEN
*  Address before (1L,JK)
                IINT = (K-KOFF)*JORB*LORB*IORB
     &               + (J-JOFF)     *LORB*IORB
     &               + (L-LOFF)          *IORB
            END IF
*
              Do i=iMin,iOff+iOrb-1
                  IJ = MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
                  KL = MAX(K,L)*(MAX(K,L)-1)/2+MIN(K,L)
                  IJKL = MAX(IJ,KL)*(MAX(IJ,KL)-1)/2+MIN(IJ,KL)
* Next line inserted by Jesper: "I don't think iInt should be the same
* for all i"
                  iInt=iInt + 1
                  Xint(iInt) = Intlst(ijkl)
              End Do
          End Do
        End Do
      End Do
*
*     Collect Exchange terms
*
      If ( IXCHNG.ne.0 ) Then
*
        iInt=0
        Do l=lOff,lOff+lOrb-1
          jMin=jOff
          If ( JLSM.ne.0 ) jMin=l
          Do j=jMin,jOff+jOrb-1
            Do k=kOff,kOff+kOrb-1
              iMin = iOff
              If(IKSM.ne.0) iMin = k
              IF(ICOUL.EQ.1)  THEN
*. Address before integral (1,j!k,l)
                  IINT = (L-LOFF)*Jorb*Korb*Iorb
     &                  + (K-KOFF)*Jorb*Iorb
     &                  + (J-JOFF)*Iorb
              ELSE IF (ICOUL.EQ.2) THEN
*  Address before (1L,JK)
                IINT = (K-KOFF)*JORB*LORB*IORB
     &               + (J-JOFF)     *LORB*IORB
     &               + (L-LOFF)          *IORB
              END IF
*
                Do i=iMin,iOff+iOrb-1
                  IL = MAX(I,L)*(MAX(I,L)-1)/2+MIN(I,L)
                  KJ = MAX(K,J)*(MAX(K,J)-1)/2+MIN(K,J)
                  ILKJ = MAX(IL,KJ)*(MAX(IL,KJ)-1)/2+MIN(IL,KJ)
                  iInt=iInt+1
                  XInt(iInt)=XInt(iInt)-Intlst(ilkj)
                End Do
            End Do
          End Do
        End Do
      End If
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(NSMOB)
      END
