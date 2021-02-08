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
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      SUBROUTINE IJIJ_CPF(JSY,HDIAG,FJI)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
      DIMENSION JSY(*),HDIAG(*),FJI(*)
      DIMENSION HCOUT(nCOP)
      JSYM(L)=JSUNP_CPF(JSY,L)
*
*
      ICOUP = 0 ! dummy initialize
      IVL   = 0 ! dummy initialize
      NSS   = 0 ! dummy initialize
*
      IADD25=IAD25S
      IAD27=0
      IF(IREF0.GT.nCOP) THEN
        WRITE(6,*)'IJIJ_CPF Error: IREF0>nCOP (See code.)'
      END IF
      CALL dDAFILE(Lu_27,2,HDIAG,IRC(1),IAD27)

      IFS=0
      TERM=0.0d0
      IVSAVE=0
      ICOUPS=0
      IOUT=0
      ICHK=0
      IADD10=IAD10(3)

300   CONTINUE
C Read a new COP buffer:
 301  CONTINUE
      CALL dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
      CALL iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
      LENGTH=ICOP1(nCOP+1)
      IF(LENGTH.EQ.0)GO TO 301
      IF(LENGTH.LT.0)GO TO 350
C A long loop over the COP buffer:
      DO 360 II=1,LENGTH
      IND=ICOP1(II)
      IF(ICHK.EQ.0) THEN
        IF(IND.NE.0)GO TO 361
        ICHK=1
        GO TO 360
      END IF

C Here, if ICHK is 1.
      ICHK=0
      INDI=IND
CPAM97      ICOUP=IAND(INDI,65535)
CPAM97      IVL=IAND(ISHFT(INDI,-16),255)
*      ICOUP=MOD(INDI,65536)
*      IVL=MOD(INDI/65536,256)
      ICOUP=IBITS(INDI,0,16)
      IVL=IBITS(INDI,16,8)
      ICHK=0
      INS=1
      IF(IVSAVE.EQ.IV0) THEN
        INS=ICOUPS
        INB=ICOUPS
      END IF

      IF(INB.GT.0) THEN
C Transfer HDIAG via buffer HCOUT, write it to unit 25:
      DO  J=INS,INB
        IOUT=IOUT+1
        HCOUT(IOUT)=HDIAG(J)
        IF(IOUT.GE.nCOP) THEN
C Write out the filled HCOUT buffer:
          IF(IFS.NE.1) THEN
            POTNUC=HCOUT(IREF0)
            IFS=1
          END IF
          DO KK=1,nCOP
            HCOUT(KK)=HCOUT(KK)-POTNUC
          END DO
          CALL dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)
          IOUT=0
        END IF
      END DO
      END IF

      IF(IVL.NE.IV0) THEN
        JJ=IRC(IVL)+ICOUP
        NSS=MUL(JSYM(JJ),LSYM)
        IF(IVL.EQ.1)INB=NVIR(NSS)
        IF(IVL.GT.1)INB=NNS(NSS)
        IF(INB.GT.0)CALL dDAFILE(Lu_27,2,HDIAG,INB,IAD27)
      END IF
      IVSAVE=IVL
      ICOUPS=ICOUP
      GO TO 360

361   CONTINUE
C Here, if ICHK.EQ.0 and IND.NE.0
CPAM97      ITYP=IAND(IND,1)
CPAM97      IJJ=IAND(ISHFT(IND,-1),2047)
*      ITYP=MOD(IND,2)
*      IJJ=MOD(IND/2,2048)
      ITYP=IBITS(IND, 0,1)
      IJJ=IBITS(IND,1,11)
      IF(ITYP.EQ.0)TERM=COP(II)*FJI(IJJ)
      IF(IVL.NE.IV0)GO TO 362

C IVL.EQ.IV0, Valence:
      INB=ICOUP
      HDIAG(INB)=HDIAG(INB)+TERM
      GO TO 360

362   IF(IVL.NE.IV1)GO TO 363
C IVL.EQ.IV1, Singles:
      INB=0
      NA1=NSYS(NSS)+1
      NA2=NSYS(NSS+1)
      IF(NA2.LT.NA1)GO TO 360
      DO NA=NA1,NA2
        INB=INB+1
        IF(ITYP.NE.0) THEN
          IIJ=IROW(LN+NA)+IJJ
          TERM=COP(II)*FJI(IIJ)
        END IF
        HDIAG(INB)=HDIAG(INB)+TERM
      END DO
      GO TO 360

363   INB=0
C Doubles:
      DO NA=1,NVIRT
        NSA=MUL(NSS,NSM(LN+NA))
        NB1=NSYS(NSA)+1
        NB2=NSYS(NSA+1)
        IF(NB2.GT.NA)NB2=NA
        IF(NB2.GE.NB1) THEN
          IIJ1=IROW(LN+NA)+IJJ
          DO NB=NB1,NB2
            INB=INB+1
            IF(ITYP.NE.0) THEN
              IIJ2=IROW(LN+NB)+IJJ
              TERM=COP(II)*(FJI(IIJ1)+FJI(IIJ2))
            END IF
            HDIAG(INB)=HDIAG(INB)+TERM
          END DO
        END IF
      END DO

360   CONTINUE
      GO TO 300

C Transfer remaining HDIAG elements to 25 via buffer HCOUT:
350   IF(INB.EQ.0)GO TO 21

      DO J=1,INB
        IOUT=IOUT+1
        HCOUT(IOUT)=HDIAG(J)
        IF(IOUT.GE.nCOP) THEN
C Write out the filled HCOUT buffer:
          IF(IFS.NE.1) THEN
            POTNUC=HCOUT(IREF0)
            IFS=1
          END IF
          DO  KK=1,nCOP
            HCOUT(KK)=HCOUT(KK)-POTNUC
          END DO
          CALL dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)
          IOUT=0
        END IF
      END DO

21    CONTINUE
C One last write of the HCOUT buffer:
      IF(IFS.NE.1) THEN
        POTNUC=HCOUT(IREF0)
        IFS=1
      END IF
      DO KK=1,IOUT
        HCOUT(KK)=HCOUT(KK)-POTNUC
      END DO
      CALL dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)
       WRITE(6,50)POTNUC
      CALL XFLUSH(6)
50    FORMAT(/,6X,'REFERENCE ENERGY',F18.8)
      RETURN
      END
