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
      SUBROUTINE RDINT2_MP2(IPRX)
C
C     SECOND ORDER TWO-ELECTRON TRANFORMATION PROGRAM. TEST SECTION
C
C     THIS SUBROUTINE READS AND CHECKS THE RESULT OF THE SECOND ORDER
C     TWO-ELECTRON TRANSFORMATION PROGRAM TRA2. IT CAN BE CALLED BY
C     TR2CTL IMMEDIATELY AFTER THE CALL TO TRA2
C
      IMPLICIT REAL*8 (A-H,O-Z)

#include "mxdim.fh"
#include "corbinf.fh"
#include "orbinf2.fh"
#include "files_mbpt2.fh"
#include "trafo.fh"
#include "WrkSpc.fh"

#include "SysDef.fh"

C
C
C     READ ADDRESS RECORD ON UNIT LUINTM
C
      IAD13=0
      LIADUT=3888
      CALL iDAFILE(LUINTM,2,IADOUT,LIADUT,IAD13)
C
C     LOOP OVER QUADRUPLES OF SYMMETRIES (NSP,NSP,NSR,NSS)
C
      ISPQRS=0
      DO 104 NSP=1,NSYM
       NBP=NBAS(NSP)
       NOP=NORB(NSP)
       NOCP=NOCC(NSP)
       DO 103 NSQ=1,NSP
        NBQ=NBAS(NSQ)
        NOQ=NORB(NSQ)
        NOCQ=NOCC(NSQ)
        NSPQ=IEOR(NSP-1,NSQ-1)+1
        DO 102 NSR=1,NSYM
         NBR=NBAS(NSR)
         NOR=NORB(NSR)
         NOCR=NOCC(NSR)
         NSPQR=IEOR(NSPQ-1,NSR-1)+1
         ISR=NSR
         DO 101 NSS=1,NSR
          NBS=NBAS(NSS)
          ISPQRS=ISPQRS+1
          IF(NSPQR.NE.NSS) GO TO 101
          NOS=NORB(NSS)
          NOCS=NOCC(NSS)
          IF(NOCP*NOCQ*NOCR*NOCS.EQ.0) GO TO 101
C
C         FIND ADDRESSES FOR THIS SYMMETRY BLOCK
C
          IADC=IADOUT(3*ISPQRS-2)
          IADX1=IADOUT(3*ISPQRS-1)
          IADX2=IADOUT(3*ISPQRS)
          WRITE(6,1000) NSP,NSQ,NSR,NSS
1000      FORMAT(/1X,'SYMMETRY BLOCK',4I4)
          IF(IADC.EQ.0) THEN
           WRITE(6,1100)
1100       FORMAT(1X,'NO COULOMB INTEGRALS FOR THIS SYMMETRY BLOCK?')
          ELSE
           WRITE(6,1200) IADC
1200       FORMAT(1X,'ADDRESS FOR COULOMB INTEGRALS',I8)
           IAD13C=IADC
          ENDIF
          IF(IADX1.EQ.0) THEN
           WRITE(6,1110)
1110       FORMAT(1X,'NO EXCHAN1 INTEGRALS FOR THIS SYMMETRY BLOCK?')
          ELSE
           WRITE(6,1210) IADX1
1210       FORMAT(1X,'ADDRESS FOR EXCHAN1 INTEGRALS',I8)
           IAD131=IADX1
          ENDIF
          IF(IADX2.EQ.0) THEN
           WRITE(6,1120)
1120       FORMAT(1X,'NO EXCHAN2 INTEGRALS FOR THIS SYMMETRY BLOCK?')
          ELSE
           WRITE(6,1220) IADX2
1220       FORMAT(1X,'ADDRESS FOR EXCHAN2 INTEGRALS',I8)
           IAD132=IADX2
          ENDIF
          LREC=NOR*NOS
          IF(NSR.EQ.NSS) LREC=(NOR**2+NOR)/2
          LRECX=NOR*NOS
          DO 10 NT=1,NOCP
           NUM=NOCQ
           IF(NSP.EQ.NSQ) NUM=NT
          DO 11 NU=1,NUM
           IF(IADC.NE.0) THEN
            Call GetMem('Tmp','ALLO','REAL',iTmp,LREC)
            CALL dDAFILE(LUINTM,2,WORK(iTmp),LREC,IAD13C)
            IF(IPRX.NE.0) LEN=LREC
            IF(IPRX.EQ.0) LEN=MIN(LREC,10)
            WRITE(6,1300) NT,NU,(WORK(I),I=iTmp,iTmp+LEN-1)
1300        FORMAT(/1X,'COULOMB INTEGRALS FOR TU PAIR',2I3
     *             /(1X,10F10.6))
            Call GetMem('Tmp','FREE','REAL',iTmp,LREC)
           ENDIF
C
C     THE LOOP ABOVE OVER T AND U RECOVERS ONE BLOCK OF INTEGRALS (AB|TU
C     FOR EACH PAIR T,U. TO PROCESS ONE BLOCK FOR ALL A AND B THE
C     FOLLOWING LOOP STRUCTURE IS USED:
C     IAB=0
C     DO 10 NA=1,NOR
C      NBM=NOS
C      IF(NSR.EQ.NSS) NBM=NA
C     DO 10 NB=1,NBM
C      IAB=IAB+1
C      WORK(IAB) NOW CONTAINS THE INTEGRAL (AB|TU)
           IF(IADX1.NE.0) THEN
            Call GetMem('Tmp','ALLO','REAL',iTmp,LRECX)
            CALL dDAFILE(LUINTM,2,WORK(iTmp),LRECX,IAD131)
            IF(IPRX.NE.0) LEN=LRECX
            IF(IPRX.EQ.0) LEN=MIN(LRECX,10)
            WRITE(6,1310) NT,NU,(WORK(I),I=iTmp,iTmp+LEN-1)
1310        FORMAT(/1X,'EXCHAN1 INTEGRALS FOR TU PAIR',2I3
     *             /(1X,10F10.6))
            Call GetMem('Tmp','FREE','REAL',iTmp,LRECX)
           ENDIF
C      THE EXCHANGE INTEGRALS OF TYPE 1 ,(AT|BU) ARE PROCESSED AS THE
C      COULOMB INTEGRALS. IF NST.NE.NSU THERE ARE ALSO EXCHANGE
C      INTEGRALS OF TYPE 2, (AU|BT). THE ORDERING IS STILL T,U AND A,B
C      BUT T IS NOW THE FOURTH INDEX AND U THE SECOND
C      EXCHANGE INTEGRALS ARE ALWAYS QUADRATIC IN A,B
           IF(IADX2.NE.0) THEN
            Call GetMem('Tmp','ALLO','REAL',iTmp,LRECX)
            CALL dDAFILE(LUINTM,2,WORK(iTmp),LRECX,IAD132)
            IF(IPRX.NE.0) LEN=LRECX
            IF(IPRX.EQ.0) LEN=MIN(LRECX,10)
            WRITE(6,1320) NT,NU,(WORK(I),I=iTmp,iTmp+LEN-1)
1320        FORMAT(/1X,'EXCHAN2 INTEGRALS FOR TU PAIR',2I3
     *             /(1X,10F10.6))
            Call GetMem('Tmp','FREE','REAL',iTmp,LRECX)
           ENDIF
11        CONTINUE
10        CONTINUE
C
C         ALL INTEGRALS FOR SYMMETRY BLOCK NSP,NSQ,NSR,NSS ARE READ
C
101      CONTINUE
102     CONTINUE
103    CONTINUE
104   CONTINUE
C
      RETURN
      END
