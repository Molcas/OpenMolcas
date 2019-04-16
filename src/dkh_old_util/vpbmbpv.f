************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1995, Bernd Artur Hess                                 *
************************************************************************
C $Id: vpbmbpv.r,v 1.4 1995/05/08 14:08:53 hess Exp $
C calculate relativistic operators
C   Bernd Artur Hess, hess@uni-bonn.de
C
C
      SUBROUTINE VPBMBPV(idbg,epsilon,EIGA,EIGB,
     *REVTA,SINVA,SINVB,NA,NB,ISIZEA,ISIZEB,VELIT,
     *AAA,AAB,RRA,RRB,PVA,VPA,ISYMA,ISYMB,BU2,G2,
     *AUX2,CMM1,BU4,CMM2,PVAA,VPAA,SCPV,SCVP)
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "real.fh"
C
      DIMENSION PVAA(ISIZEA),VPAA(ISIZEA)
C
      DIMENSION EIGA(NA,NA),EIGB(NB,NB),SINVA(NA,NA),
     *SINVB(NB,NB),REVTA(NA,NA),
     *AAA(NA),RRA(NA),AAB(NB),RRB(NB),
     *PVA(NA,NB),VPA(NA,NB),
     *BU2(NA,NB),G2(NA,NB),
     *CMM1(NA,NB),BU4(NA,NB)
C
      DIMENSION AUX2(NA,NB),CMM2(NA,NB),
     *          SCPV(NB,NA),SCVP(NB,NA)
C
*     Call qEnter('vpbmbpv')
C
C
      IF(iSyma.eq.iSymb) Then
        IJ=0
        DO I=1,NA
        DO J=1,I
        IJ=IJ+1
        PVA(I,J)=PVAA(IJ)
        PVA(J,I)=-VPAA(IJ)
        VPA(I,J)=VPAA(IJ)
        VPA(J,I)=-PVAA(IJ)
        ENDDO
        ENDDO
      ENDIF
C
      IF(iSyma.lt.iSymb) THEN
      DO I=1,NA
      DO J=1,NB
      CMM1(I,J) = SCPV(J,I)
      CMM2(I,J) = SCVP(J,I)
      ENDDO
      ENDDO
      DO I=1,NA
      DO J=1,NB
      PVA(I,J) = CMM1(I,J)
      VPA(I,J) = CMM2(I,J)
      ENDDO
      ENDDO
      CALL DCOPY_(NA*NB,[ZERO],0,CMM1,1)
      CALL DCOPY_(NA*NB,[ZERO],0,CMM2,1)
C
      ENDIF
C
C
C    TRANSFORM pV TO T-BASIS
C
C
      CALL TrSmrN(PVA,SINVA,SINVB,G2,NA,NB,AUX2,CMM1)
      CALL TrSmrN(G2,EIGA,EIGB,BU2,NA,NB,AUX2,CMM1)
*
*
      call dcopy_(na*nb,[Zero],0,G2,1)
      call dcopy_(na*nb,[Zero],0,Aux2,1)
*
C    TRANSFORM Vp TO T-BASIS
*
      CALL TrSmrN(VPA,SINVA,SINVB,G2,NA,NB,AUX2,CMM1)
      CALL TrSmrN(G2,EIGA,EIGB,BU4,NA,NB,AUX2,CMM1)
*
*
C
      call dcopy_(na*nb,[Zero],0,G2,1)
      call dcopy_(na*nb,[Zero],0,Aux2,1)
      call dcopy_(na*nb,[Zero],0,Cmm1,1)
C
C     Multiply
C
C
      DO I=1,NA
      DO J=1,NB
      G2(I,J)=BU4(I,J)*RRB(J)
      CMM1(I,J)=-BU2(I,J)*RRA(I)
      Enddo
      Enddo
C
*
*
*     write(*,*)
*     write(*,*) 'pV part ofcommutator MATRIX'
*     write(*,*)  G2
*     write(*,*)
*
*     write(*,*)
*     write(*,*) '-Vp part commutator MATRIX'
*     write(*,*)  CMM1
*     write(*,*)
*
*
C     Caculate the commutator <iSymA|[V,pb]|iSymB> and put into CMM1
C
      CAll AddMar(NA*NB,G2,CMM1)
C
*
C
C     MultiplY BY A MATRIX
C
      DO I=1,NA
      DO J=1,NB
      CMM1(I,J)=CMM1(I,J)*AAA(I)*AAB(J)
      ENDDO
      ENDDO
C
C
C     MultiplY BY REVTA  MATRIX
C
C
      DO I=1,NA
      DO J=1,NB
      DO K=1,NA
      CMM2(I,J)=CMM2(I,J)+REVTA(I,K)*CMM1(K,J)
      ENDDO
      ENDDO
      ENDDO
C
      call dcopy_(na*nb,[Zero],0,Cmm1,1)
C
*     write(*,*)
*     write(*,*) 'CMM2  MATRX FINAL'
*     write(*,*)  CMM2
*     write(*,*)
C
C
*     write(*,*) 'END OF VPBMBPV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
C
*     Call qExit('vpbmbpv')
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(idbg)
        CALL Unused_real(epsilon)
        CALL Unused_integer(ISIZEB)
        CALL Unused_real(VELIT)
      END IF
      END
      SUBROUTINE TrSmrN(A,BA,BB,C,NA,NB,H,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NA,NB),BA(NA,NA),BB(NB,NB),
     &C(NA,NB),H(NA,NB),W(NA,NB)
C
C     TRANSFORM SYMMETRIC MATRIX A BY UNITARY TRANSFORMATION
C     IN B. RESULT IS IN C
C
C
      CALL DZero(C,NA*NB)
      CALL DZero(H,NA*NB)
      call dcopy_(NA*NB,A,1,W,1)
*
C
      DO 2 I=1,NA
      DO 3 L=1,NB
      DO 4 K=1,NA
      H(I,L)=BA(K,I)*W(K,L)+H(I,L)
4     CONTINUE
3     CONTINUE
2     CONTINUE
*
*
      DO 21 I=1,NA
      DO 22 J=1,NB
      DO 23 L=1,NB
      C(I,J)=H(I,L)*BB(L,J)+C(I,J)
23    CONTINUE
22    CONTINUE
21    CONTINUE
*
*
      RETURN
      END
