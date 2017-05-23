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
      SUBROUTINE IJIJ(INTSYM,HDIAG,FC,FIJIJ)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION INTSYM(*),HDIAG(*),FC(*),
     *          FIJIJ(*)
      DIMENSION HCOUT(nCOP)
*
      JSYM(L)=JSUNP(INTSYM,L)
*------
* POW: Unnecessary but warning stopping initializations
      inb=-1234567
*------
      IADD25=IAD25S
      IAD27=0
      IREF0=1
      CALL dDAFILE(Lu_27,2,HDIAG,IRC(1),IAD27)
*
*     WRITE(6,*) ' Hdiag'
*     WRITE(6,*) ( Hdiag(i),i=1,IRC(1) )
*
      IFS=0
      IVL=0
      IVSAVE=0
      ICOUPS=0
      ICOUP =0
      NSS=1
      IOUT=0
      ICHK=0
      IADD10=IAD10(3)
      TERM=0.0D00
300   CONTINUE
C READ A COP BUFFER:
      CALL dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
      CALL iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
      LENGTH=ICOP1(nCOP+1)
      IF(LENGTH.EQ.0)GO TO 300
      IF(LENGTH.LT.0)GO TO 350
C LOOP OVER THE COP BUFFER:
      DO 360 II=1,LENGTH
      IND=ICOP1(II)
      IF(ICHK.NE.0)GO TO 460
      IF(IND.NE.0)GO TO 361
      ICHK=1
      GO TO 360
460   ICHK=0
      INDI=IND
*      ICOUP=MOD(INDI,2**16)
*      IVL=MOD(INDI/2**16,2**8)
      ICOUP=IBITS(INDI, 0,16)
      IVL=IBITS(INDI,16,8)
      ICHK=0
      INS=1
      IF(IVSAVE.EQ.IVVER) THEN
        INS=ICOUPS
        INB=ICOUPS
      END IF
      IF(INB.NE.0) THEN
        DO 10 J=INS,INB
          IOUT=IOUT+1
          HCOUT(IOUT)=HDIAG(J)
          IF(IOUT.LT.nCOP)GO TO 10
          IF(IFS.EQ.0)THEN
            POTNUC=HCOUT(IREF0)
            IFS=1
          END IF
          DO 8831 KK=1,nCOP
            HCOUT(KK)=HCOUT(KK)-POTNUC
8831      CONTINUE
          CALL dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)
          IOUT=0
10      CONTINUE
      END IF
      IF(IVL.NE.IVVER) THEN
        JJ=IRC(IVL)+ICOUP
        NSS=MUL(JSYM(JJ),LSYM)
        IF(IVL.EQ.IDVER)THEN
          INB=NVIR(NSS)
        ELSE
          INB=NVPAIR(NSS)
        END IF
        IF(INB.GT.0)CALL dDAFILE(Lu_27,2,HDIAG,INB,IAD27)
      END IF
      IVSAVE=IVL
      ICOUPS=ICOUP
      GO TO 360
361   CONTINUE
*      ITYP=MOD(IND,2)
*      IJJ=MOD(IND/2,2**11)
      ITYP=IBITS(IND,0,1)
      IJJ=IBITS(IND,1,11)
      IF(ITYP.EQ.0)TERM=COP(II)*FIJIJ(IJJ)
      IF(IVL.EQ.IVVER) THEN
        INB=ICOUP
        HDIAG(INB)=HDIAG(INB)+TERM
      ELSE IF(IVL.EQ.IDVER) THEN
        INB=0
        NA1=NVIRP(NSS)+1
        NA2=NVIRP(NSS)+NVIR(NSS)
        IF(NA2.LT.NA1)GO TO 360
        DO 365 NA=NA1,NA2
          INB=INB+1
          IF(ITYP.EQ.1) THEN
            IIJ=IROW(LN+NA)+IJJ
            TERM=COP(II)*FIJIJ(IIJ)
          END IF
          HDIAG(INB)=HDIAG(INB)+TERM
365     CONTINUE
      ELSE
        INB=0
        DO 375 NA=1,NVIRT
          NSA=MUL(NSS,NSM(LN+NA))
          NB1=NVIRP(NSA)+1
          NB2=NVIRP(NSA)+NVIR(NSA)
          IF(NB2.GT.NA)NB2=NA
          IF(NB2.LT.NB1)GO TO 375
          IIJ1=IROW(LN+NA)+IJJ
          DO 376 NB=NB1,NB2
            INB=INB+1
            IF(ITYP.EQ.1) THEN
              IIJ2=IROW(LN+NB)+IJJ
              TERM=COP(II)*(FIJIJ(IIJ1)+FIJIJ(IIJ2))
            END IF
            HDIAG(INB)=HDIAG(INB)+TERM
376       CONTINUE
375     CONTINUE
      END IF
360   CONTINUE
      GO TO 300
C EMPTY LAST BUFFER
350   CONTINUE
      DO 20 J=1,INB
        IOUT=IOUT+1
        HCOUT(IOUT)=HDIAG(J)
        IF(IOUT.LT.nCOP)GO TO 20
        IF(IFS.EQ.0)THEN
          POTNUC=HCOUT(IREF0)
          IFS=1
        END IF
        DO 8830 KK=1,nCOP
          HCOUT(KK)=HCOUT(KK)-POTNUC
8830    CONTINUE
        CALL dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)
        IOUT=0
20    CONTINUE
      IF(IFS.EQ.0)THEN
        POTNUC=HCOUT(IREF0)
        IFS=1
      END IF
      DO 8829 KK=1,IOUT
        HCOUT(KK)=HCOUT(KK)-POTNUC
8829  CONTINUE
      CALL dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_real_array(FC)
      END
