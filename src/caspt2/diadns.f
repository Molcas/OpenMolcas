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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*

      SUBROUTINE DIADNS(ISYM,ICASE,VEC1,VEC2,DPT2,LIST)

      IMPLICIT REAL*8 (A-H,O-Z)


#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "sigma.fh"

      DIMENSION VEC1(*),VEC2(*)
      DIMENSION DPT2(*)
      DIMENSION LIST(*)
      DIMENSION IOFDIJ(8),IOFDAB(8)
      DIMENSION IOFCD(8,8)

C Compute diagonal-block contribs to a trans density matrix.
C Vector blocks are in spectral resolution basis (ON).
C Each square matrix block of density matrix elements is
C computed and stored in full, even if VEC1=VEC2.
C Present implementation does not compute active-active
C contributions. This should be added in a separate routine,
C since it requires transformation to standard (Non-ON) basis.

      NIN=NINDEP(ISYM,ICASE)
      IF(NIN.EQ.0) RETURN
      NIS=NISUP(ISYM,ICASE)
      IF(NIS.EQ.0) RETURN
      CALL QENTER('DIADNS')
      NAS=NASUP(ISYM,ICASE)
      NVEC=NIN*NIS
      SQR2=SQRT(2.0D00)

      IFTEST=0
C Set up various offset arrays:
      IDIJ=0
      DO IS=1,NSYM
        NI=NISH(IS)
        NA=NASH(IS)
        NO=NORB(IS)
        IDTU=IDIJ+NO*NI+NI
        IDAB=IDTU+NO*NA+NA
        IOFDIJ(IS)=IDIJ
        IOFDAB(IS)=IDAB
        IDIJ=IDIJ+NO*NO
      END DO
      DO IS=1,NSYM
        ICD=0
        ICEP=0
        ICEM=0
        ICGP=0
        ICGM=0
        DO JS=1,NSYM
          IJS=MUL(IS,JS)
          IOFCD(IS,JS)=ICD
          ICD =ICD +NSSH(JS)*NISH(IJS)
          ICEP=ICEP+NSSH(JS)*NIGEJ(IJS)
          ICEM=ICEM+NSSH(JS)*NIGTJ(IJS)
          ICGP=ICGP+NISH(JS)*NAGEB(IJS)
          ICGM=ICGM+NISH(JS)*NAGTB(IJS)
        END DO
      END DO

C Core contribution:
      OVL=DDOT_(NVEC,VEC1,1,VEC2,1)
      DO IS=1,NSYM
        NI=NISH(IS)
        NO=NORB(IS)
        IDII=IOFDIJ(IS)+1
        DO III=1,NI
          DPT2(IDII)=DPT2(IDII)+2.0D0*OVL
          IDII=IDII+NO+1
        END DO
      END DO
*
      LLST1 = 0 ! dummy initialize
      NLST1 = 0 ! dummy initialize
*
      GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13) ICASE
      CALL QEXIT('DIADNS')
      RETURN

C -----------------------------------------------
   1  CONTINUE
C Case A
      NI=NISH(ISYM)
      NO=NORB(ISYM)
      DO II=1,NI
        IV2=1+NIN*(II-1)
        DO IJ=1,NI
          IDIJ=IOFDIJ(ISYM)+II+NO*(IJ-1)
          IV1=1+NIN*(IJ-1)
          DPT2(IDIJ)=DPT2(IDIJ)-
     &           DDOT_(NIN,VEC1(IV1),1,VEC2(IV2),1)
        END DO
      END DO
      GOTO 100
C -----------------------------------------------
   2  CONTINUE
C Case BP
   3  CONTINUE
C Case BM
C Unfold VEC1 and VEC2 into X1(MU,K,I), X2(MU,K,I):
      NX=NIN*NIMX**2
      CALL GETMEM('DIA_X1','ALLO','REAL',LX1,NX)
      CALL GETMEM('DIA_X2','ALLO','REAL',LX2,NX)
      DO ISYMK=1,NSYM
       ISYMI=MUL(ISYMK,ISYM)
       NK=NISH(ISYMK)
       NI=NISH(ISYMI)
       NKI=NK*NI
       IF(NKI.EQ.0) GOTO 903
       CALL DCOPY_(NIN*NKI,[0.0D0],0,WORK(LX1),1)
       CALL DCOPY_(NIN*NKI,[0.0D0],0,WORK(LX2),1)
       IF(ICASE.EQ.2) THEN
         LLST1=LLIST(ISYMK,ISYM,14)
         NLST1=NLIST(ISYMK,ISYM,14)
         VAL1(1)= 1.0D00
         VAL1(2)= SQR2
       ELSE IF(ICASE.EQ.3) THEN
         LLST1=LLIST(ISYMK,ISYM,15)
         NLST1=NLIST(ISYMK,ISYM,15)
         VAL1(1)= 1.0D00
         VAL1(2)=-1.0D00
       END IF
       IF(NLST1.EQ.0) GOTO 903
       INCX1=1
       INCX2=NIN
       INCX3=NIN*NK
       INCY1=1
       INCY2=NIN
       LEN1=NIN
       CALL MLTUNF(LIST(LLST1),WORK(LX1),VEC1)
       CALL MLTUNF(LIST(LLST1),WORK(LX2),VEC2)
C D(I,J) := Add contraction -X2(MU,K,I)*X1(MU,K,J):
       IDIJ=1+IOFDIJ(ISYMI)
       NO=NORB(ISYMI)
       CALL DGEMM_('T','N',NI,NI,NIN*NK,-1.0D00,
     &            WORK(LX2),NIN*NK,WORK(LX1),NIN*NK,
     &            1.0D00,DPT2(IDIJ),NO)
 903   CONTINUE
      END DO
      CALL GETMEM('DIA_X1','FREE','REAL',LX1,NX)
      CALL GETMEM('DIA_X2','FREE','REAL',LX2,NX)
      GOTO 100
C -----------------------------------------------
   4  CONTINUE
C Case C
      NS=NSSH(ISYM)
      NO=NORB(ISYM)
      DO IA=1,NS
        IV1=1+NIN*(IA-1)
        DO IB=1,NS
          IDAB=IOFDAB(ISYM)+IA+NO*(IB-1)
          IV2=1+NIN*(IB-1)
          DPT2(IDAB)=DPT2(IDAB)+
     &          DDOT_(NIN,VEC1(IV1),1,VEC2(IV2),1)
        END DO
      END DO
      GOTO 100
C -----------------------------------------------
   5  CONTINUE
C Case D
      DO ISYMA=1,NSYM
       NS=NSSH(ISYMA)
       NOA=NORB(ISYMA)
       ISYMI=MUL(ISYMA,ISYM)
       NI=NISH(ISYMI)
       NOI=NORB(ISYMI)
       IV=1+NIN*IOFCD(ISYM,ISYMA)
       INCA=NIN*NI
       DO II=1,NI
         IV2=IV+NIN*(II-1)
         DO IJ=1,NI
           IDIJ=IOFDIJ(ISYMI)+II+NOI*(IJ-1)
           SUM=DPT2(IDIJ)
           IV1=IV+NIN*(IJ-1)
           DO IA=1,NS
             IV11=IV1+INCA*(IA-1)
             IV22=IV2+INCA*(IA-1)
             SUM=SUM-DDOT_(NIN,VEC1(IV11),1,VEC2(IV22),1)
           END DO
           DPT2(IDIJ)=SUM
         END DO
       END DO
       DO IA=1,NS
         IV1=IV+INCA*(IA-1)
         DO IB=1,NS
           IDAB=IOFDAB(ISYMA)+IA+NOA*(IB-1)
           IV2=IV+INCA*(IB-1)
             DPT2(IDAB)=DPT2(IDAB)+
     &             DDOT_(INCA,VEC1(IV1),1,VEC2(IV2),1)
         END DO
       END DO
      END DO

      GOTO 100
C -----------------------------------------------
   6  CONTINUE
C Case EP
   7  CONTINUE
C Case EM
      NX=NIN*NSMX*NIMX**2
      CALL GETMEM('DIA_X1','ALLO','REAL',LX1,NX)
      CALL GETMEM('DIA_X2','ALLO','REAL',LX2,NX)
      IYOFF=0
      DO ISYMA=1,NSYM
       ISYMKI=MUL(ISYMA,ISYM)
       NA=NSSH(ISYMA)
       NOA=NORB(ISYMA)
       IF(ICASE.EQ.6) NKIY=NIGEJ(ISYMKI)
       IF(ICASE.EQ.7) NKIY=NIGTJ(ISYMKI)
       IY=1+IYOFF
C First, contributions to DIJ.
C Unfold VEC1 and VEC2 into X1(MU,A;K,I), X2(MU,A;K,I):
       DO ISYMK=1,NSYM
        ISYMI=MUL(ISYMK,ISYMKI)
        NK=NISH(ISYMK)
        NI=NISH(ISYMI)
        NAKI=NA*NK*NI
        IF(NAKI.EQ.0) GOTO 907
        CALL DCOPY_(NIN*NAKI,[0.0D0],0,WORK(LX1),1)
        CALL DCOPY_(NIN*NAKI,[0.0D0],0,WORK(LX2),1)
        IF(ICASE.EQ.6) THEN
          LLST1=LLIST(ISYMK,ISYMKI,14)
          NLST1=NLIST(ISYMK,ISYMKI,14)
          VAL1(1)= 1.0D00
          VAL1(2)= SQR2
        ELSE IF(ICASE.EQ.7) THEN
          LLST1=LLIST(ISYMK,ISYMKI,15)
          NLST1=NLIST(ISYMK,ISYMKI,15)
          VAL1(1)= 1.0D00
          VAL1(2)=-1.0D00
        END IF
        IF(NLST1.EQ.0) GOTO 907
        INCX1=1
        INCX2=NIN*NA
        INCX3=NIN*NA*NK
        INCY1=1
        INCY2=NIN*NA
        LEN1=NIN*NA
        CALL MLTUNF(LIST(LLST1),WORK(LX1),VEC1(IY))
        CALL MLTUNF(LIST(LLST1),WORK(LX2),VEC2(IY))
C  D(I,J) := Add contraction -X2(MU,A,K,I)*X1(MU,A,K,J):
        IDIJ=1+IOFDIJ(ISYMI)
        NOI=NORB(ISYMI)
        CALL DGEMM_('T','N',NI,NI,NIN*NA*NK,-1.0D00,
     &             WORK(LX2),NIN*NA*NK,WORK(LX1),NIN*NA*NK,
     &             1.0D00,DPT2(IDIJ),NOI)
 907    CONTINUE
       END DO
C Second, contributions to DAB.
       IF(NKIY.GT.0) THEN
        DO IA=1,NA
         DO IB=1,NA
          IDAB=IOFDAB(ISYMA)+IA+NOA*(IB-1)
          SUM=DPT2(IDAB)
          DO MU=1,NIN
            IY1=IYOFF+MU+NIN*(IA-1)
            IY2=IYOFF+MU+NIN*(IB-1)
            SUM=SUM+DDOT_(NKIY,VEC1(IY1),INCY2,VEC2(IY2),INCY2)
          END DO
          DPT2(IDAB)=SUM
         END DO
        END DO
       END IF
       IYOFF=IYOFF+NIN*NA*NKIY
      END DO
      CALL GETMEM('DIA_X1','FREE','REAL',LX1,NX)
      CALL GETMEM('DIA_X2','FREE','REAL',LX2,NX)
      GOTO 100
C -----------------------------------------------
   8  CONTINUE
C Case FP
   9  CONTINUE
C Case FM
C Unfold VEC1 and VEC2 into X1(MU,C,A), X2(MU,C,B):
      NX=NIN*NSMX**2
      CALL GETMEM('DIA_X1','ALLO','REAL',LX1,NX)
      CALL GETMEM('DIA_X2','ALLO','REAL',LX2,NX)
      DO ISYMC=1,NSYM
       ISYMA=MUL(ISYMC,ISYM)
       NC=NSSH(ISYMC)
       NA=NSSH(ISYMA)
       NCA=NC*NA
       IF(NCA.EQ.0) GOTO 909
       CALL DCOPY_(NIN*NCA,[0.0D0],0,WORK(LX1),1)
       CALL DCOPY_(NIN*NCA,[0.0D0],0,WORK(LX2),1)
       IF(ICASE.EQ.8) THEN
         LLST1=LLIST(ISYMC,ISYM,16)
         NLST1=NLIST(ISYMC,ISYM,16)
         VAL1(1)= 1.0D00
         VAL1(2)= SQR2
       ELSE IF(ICASE.EQ.9) THEN
         LLST1=LLIST(ISYMC,ISYM,17)
         NLST1=NLIST(ISYMC,ISYM,17)
         VAL1(1)= 1.0D00
         VAL1(2)=-1.0D00
       END IF
       IF(NLST1.EQ.0) GOTO 909
       INCX1=1
       INCX2=NIN
       INCX3=NIN*NC
       INCY1=1
       INCY2=NIN
       LEN1=NIN
       CALL MLTUNF(LIST(LLST1),WORK(LX1),VEC1)
       CALL MLTUNF(LIST(LLST1),WORK(LX2),VEC2)
C D(A,B) := Add contraction  X1(MU,C,A)*X2(MU,C,B):
       IDAB=1+IOFDAB(ISYMA)
       NOA=NORB(ISYMA)
       CALL DGEMM_('T','N',NA,NA,NIN*NC,+1.0D00,
     &            WORK(LX1),NIN*NC,WORK(LX2),NIN*NC,
     &            1.0D00,DPT2(IDAB),NOA)
 909   CONTINUE
      END DO
      CALL GETMEM('DIA_X1','FREE','REAL',LX1,NX)
      CALL GETMEM('DIA_X2','FREE','REAL',LX2,NX)
      GOTO 100
C -----------------------------------------------
  10  CONTINUE
C Case GP
  11  CONTINUE
C Case GM
      NX=NIN*NIMX*NSMX**2
      CALL GETMEM('DIA_X1','ALLO','REAL',LX1,NX)
      CALL GETMEM('DIA_X2','ALLO','REAL',LX2,NX)
      IYOFF=0
      DO ISYMI=1,NSYM
       ISYMCA=MUL(ISYMI,ISYM)
       NI=NISH(ISYMI)
       NOI=NORB(ISYMI)
       IF(ICASE.EQ.10) NCAY=NAGEB(ISYMCA)
       IF(ICASE.EQ.11) NCAY=NAGTB(ISYMCA)
       IY=1+IYOFF
C First, contributions to DAB.
C Unfold VEC1 and VEC2 into X1(MU,I;C,A), X2(MU,I;C,A):
       DO ISYMC=1,NSYM
        ISYMA=MUL(ISYMC,ISYMCA)
        NC=NSSH(ISYMC)
        NA=NSSH(ISYMA)
        NICA=NI*NC*NA
        IF(NICA.EQ.0) GOTO 911
        CALL DCOPY_(NIN*NICA,[0.0D0],0,WORK(LX1),1)
        CALL DCOPY_(NIN*NICA,[0.0D0],0,WORK(LX2),1)
        IF(ICASE.EQ.10) THEN
          LLST1=LLIST(ISYMC,ISYMCA,16)
          NLST1=NLIST(ISYMC,ISYMCA,16)
          VAL1(1)= 1.0D00
          VAL1(2)= SQR2
        ELSE IF(ICASE.EQ.11) THEN
          LLST1=LLIST(ISYMC,ISYMCA,17)
          NLST1=NLIST(ISYMC,ISYMCA,17)
          VAL1(1)= 1.0D00
          VAL1(2)=-1.0D00
        END IF
        IF(NLST1.EQ.0) GOTO 911
        INCX1=1
        INCX2=NIN*NI
        INCX3=NIN*NI*NC
        INCY1=1
        INCY2=NIN*NI
        LEN1=NIN*NI
        CALL MLTUNF(LIST(LLST1),WORK(LX1),VEC1(IY))
        CALL MLTUNF(LIST(LLST1),WORK(LX2),VEC2(IY))
C  D(A,B) := Add contraction +X1(MU,I,C,A)*X2(MU,I,C,B):
        IDAB=1+IOFDAB(ISYMA)
        NOA=NORB(ISYMA)
        CALL DGEMM_('T','N',NA,NA,NIN*NI*NC,+1.0D00,
     &             WORK(LX1),NIN*NI*NC,WORK(LX2),NIN*NI*NC,
     &             1.0D00,DPT2(IDAB),NOA)
 911    CONTINUE
       END DO
C Second, contributions to DIJ.
       IF(NCAY.GT.0) THEN
        DO II=1,NI
         DO IJ=1,NI
          IDIJ=IOFDIJ(ISYMI)+II+NOI*(IJ-1)
          SUM=DPT2(IDIJ)
          DO MU=1,NIN
            IY1=IYOFF+MU+NIN*(IJ-1)
            IY2=IYOFF+MU+NIN*(II-1)
            SUM=SUM-DDOT_(NCAY,VEC1(IY1),INCY2,VEC2(IY2),INCY2)
          END DO
          DPT2(IDIJ)=SUM
         END DO
        END DO
       END IF
       IYOFF=IYOFF+NIN*NI*NCAY
      END DO
      CALL GETMEM('DIA_X1','FREE','REAL',LX1,NX)
      CALL GETMEM('DIA_X2','FREE','REAL',LX2,NX)
      GOTO 100
C -----------------------------------------------
  12  CONTINUE
C Case HP
  13  CONTINUE
C Case HM
C Unfold VEC1 and VEC2 into X1(MU,K,I), X2(MU,K,I):
      NX=NAS*NIMX**2
      CALL GETMEM('DIA_X1','ALLO','REAL',LX1,NX)
      CALL GETMEM('DIA_X2','ALLO','REAL',LX2,NX)
      DO ISYMK=1,NSYM
       ISYMI=MUL(ISYMK,ISYM)
       NK=NISH(ISYMK)
       NI=NISH(ISYMI)
       NKI=NK*NI
       IF(NKI.EQ.0) GOTO 813
       CALL DCOPY_(NAS*NKI,[0.0D0],0,WORK(LX1),1)
       CALL DCOPY_(NAS*NKI,[0.0D0],0,WORK(LX2),1)
       IF(ICASE.EQ.12) THEN
         LLST1=LLIST(ISYMK,ISYM,14)
         NLST1=NLIST(ISYMK,ISYM,14)
         VAL1(1)= 1.0D00
         VAL1(2)= SQR2
       ELSE IF(ICASE.EQ.13) THEN
         LLST1=LLIST(ISYMK,ISYM,15)
         NLST1=NLIST(ISYMK,ISYM,15)
         VAL1(1)= 1.0D00
         VAL1(2)=-1.0D00
       END IF
       IF(NLST1.EQ.0) GOTO 813
       INCX1=1
       INCX2=NAS
       INCX3=NAS*NK
       INCY1=1
       INCY2=NAS
       LEN1=NAS
       CALL MLTUNF(LIST(LLST1),WORK(LX1),VEC1)
       CALL MLTUNF(LIST(LLST1),WORK(LX2),VEC2)
C D(I,J) := Add contraction -X2(MU,K,I)*X1(MU,K,J):
       IDIJ=1+IOFDIJ(ISYMI)
       NOI=NORB(ISYMI)
       CALL DGEMM_('T','N',NI,NI,NAS*NK,-1.0D00,
     &            WORK(LX2),NAS*NK,WORK(LX1),NAS*NK,
     &            1.0D00,DPT2(IDIJ),NOI)
 813   CONTINUE
      END DO
      CALL GETMEM('DIA_X1','FREE','REAL',LX1,NX)
      CALL GETMEM('DIA_X2','FREE','REAL',LX2,NX)
C Unfold VEC1 and VEC2 into X1(A,C,IJ), X2(A,C,IJ):
      NX=NIS*NSMX**2
      CALL GETMEM('DIA_X1','ALLO','REAL',LX1,NX)
      CALL GETMEM('DIA_X2','ALLO','REAL',LX2,NX)
      DO ISYMC=1,NSYM
       ISYMA=MUL(ISYMC,ISYM)
       NC=NSSH(ISYMC)
       NA=NSSH(ISYMA)
       NCA=NC*NA
       IF(NCA.EQ.0) GOTO 913
       CALL DCOPY_(NIS*NCA,[0.0D0],0,WORK(LX1),1)
       CALL DCOPY_(NIS*NCA,[0.0D0],0,WORK(LX2),1)
       IF(ICASE.EQ.12) THEN
         LLST1=LLIST(ISYMA,ISYM,16)
         NLST1=NLIST(ISYMA,ISYM,16)
         VAL1(1)= 1.0D00
         VAL1(2)= SQR2
       ELSE IF(ICASE.EQ.13) THEN
         LLST1=LLIST(ISYMA,ISYM,17)
         NLST1=NLIST(ISYMA,ISYM,17)
         VAL1(1)= 1.0D00
         VAL1(2)=-1.0D00
       END IF
       IF(NLST1.EQ.0) GOTO 913
       INCX1=NCA
       INCX2=1
       INCX3=NA
       INCY1=NAS
       INCY2=1
       LEN1=NIS
       CALL MLTUNF(LIST(LLST1),WORK(LX1),VEC1)
       CALL MLTUNF(LIST(LLST1),WORK(LX2),VEC2)
C D(A,B) := Add contraction  X1(A,C,IJ)*X2(B,C,IJ):
       IDAB=1+IOFDAB(ISYMA)
       NOA=NORB(ISYMA)
       CALL DGEMM_('N','T',NA,NA,NIS*NC,+1.0D00,
     &            WORK(LX1),NA,WORK(LX2),NA,
     &            1.0D00,DPT2(IDAB),NOA)
 913   CONTINUE
      END DO
      CALL GETMEM('DIA_X1','FREE','REAL',LX1,NX)
      CALL GETMEM('DIA_X2','FREE','REAL',LX2,NX)
      GOTO 100
C -----------------------------------------------
 100  CONTINUE

      CALL QEXIT('DIADNS')
      RETURN
      END
