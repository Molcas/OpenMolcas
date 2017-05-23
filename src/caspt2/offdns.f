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
      SUBROUTINE OFFDNS(ISYM1,ICASE1,ISYM2,ICASE2,X1,X2,DPT2,Y,LIST)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
      DIMENSION X1(*),X2(*),Y(*)
      DIMENSION DPT2(*)
      DIMENSION LIST(*)
      DIMENSION IOFDIT(8),IOFDIA(8),IOFDTA(8)
      DIMENSION IOFCD(8,8),IOFCEP(8,8),IOFCEM(8,8),IOFCGP(8,8),
     &          IOFCGM(8,8)
#include "sigma.fh"
      DIMENSION IFCOUP(13,13)
      DATA IFCOUP / 0, 1, 2, 0, 3, 4, 5, 0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 0, 8, 0, 0, 9,10,11,12, 0, 0,
     &              0, 0, 0, 0, 0,13,14, 0, 0,15,16,23,24,
     &              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,17, 0,
     &              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,18,
     &              0, 0, 0, 0, 0, 0, 0, 0, 0,19, 0, 0, 0,
     &              0, 0, 0, 0, 0, 0, 0, 0, 0, 0,20, 0, 0,
     &              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,21, 0,
     &              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,22,
     &              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
C Compute off-diagonal contributions to a trans density matrix.
C Sub-diagonal blocks only. If a density matrix is required,
C i.e. both wave functions equal, use symmetry. Else, call
C again with wave functions interchanged.

C Various constants:
      SQR2=SQRT(2.0D00)
      SQR3=SQRT(3.0D00)
      SQR6=SQRT(6.0D00)
      SQRI2=1.0D00/SQR2
      SQRI6=1.0D00/SQR6
      SQR32=SQR3*SQRI2

*
      NA = 0 ! dummy initialize
*
      IFTEST=0
      IMLTOP=2
      KOD=IFCOUP(ICASE2,ICASE1)
      IF(KOD.EQ.0) RETURN

CPAM      IF(KOD.EQ.23 .OR. KOD.EQ.24) IFTEST=1
CPAM      if(iftest.eq.1) then
CPAM      WRITE(*,*)' OFFDNS ICASE1=',ICASE1
CPAM      WRITE(*,*)'        ICASE2=',ICASE2
CPAM      WRITE(*,*)'        IFTEST=',IFTEST
CPAM      end if

      ISYM12=MUL(ISYM1,ISYM2)
      NAS1=NASUP(ISYM1,ICASE1)
      NIS1=NISUP(ISYM1,ICASE1)
      NAS2=NASUP(ISYM2,ICASE2)
C Set up various offset arrays:
      IDOFF=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NO=NORB(ISYM)
        IOFDIT(ISYM)=IDOFF+NO*NI
        IOFDIA(ISYM)=IOFDIT(ISYM)+NO*NA
        IOFDTA(ISYM)=IOFDIA(ISYM)+NI
        IDOFF=IDOFF+NO**2
        ICD=0
        ICEP=0
        ICEM=0
        ICGP=0
        ICGM=0
        DO JSYM=1,NSYM
          IJSYM=MUL(ISYM,JSYM)
          IOFCD(ISYM,JSYM)=ICD
          IOFCEP(ISYM,JSYM)=ICEP
          IOFCEM(ISYM,JSYM)=ICEM
          IOFCGP(ISYM,JSYM)=ICGP
          IOFCGM(ISYM,JSYM)=ICGM
          ICD =ICD +NSSH(JSYM)*NISH(IJSYM)
          ICEP=ICEP+NSSH(JSYM)*NIGEJ(IJSYM)
          ICEM=ICEM+NSSH(JSYM)*NIGTJ(IJSYM)
          ICGP=ICGP+NISH(JSYM)*NAGEB(IJSYM)
          ICGM=ICGM+NISH(JSYM)*NAGTB(IJSYM)
        END DO
      END DO

      GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     &      16,17,18,19,20,21,22,23,24) KOD
      RETURN

C  -----------------------------------------------
   1  CONTINUE
C  A&BP One-el
      LLST1=LLIST(ISYM1,ISYM2,12)
      NLST1=NLIST(ISYM1,ISYM2,12)
      IF(NLST1.EQ.0) GOTO 91
      VAL1(1)=1.0D00
      VAL1(2)=2.0D00
      LLST2=LLIST(ISYM1,ISYM2,14)
      NLST2=NLIST(ISYM1,ISYM2,14)
      IF(NLST2.EQ.0) GOTO 91
      VAL2(1)=1.0D00
      VAL2(2)=SQR2
      IXTI=1
      INCX1=1
      INCX2=NASH(ISYM1)
      IDIT=1+IOFDIT(ISYM12)
      INCF1=NORB(ISYM12)
      INCF2=1
      IY=1
      INCY1=1
      INCY2=NTGEU(ISYM2)
      CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &            X1(IXTI),DPT2(IDIT),Y(IY))
  91  CONTINUE

C  A&BP Two-el
      LLST1=LLIST(ISYM1,ISYM2, 3)
      NLST1=NLIST(ISYM1,ISYM2, 3)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)=-1.0D00
      VAL1(2)=-2.0D00
      LLST2=LLIST(ISYM1,ISYM2,14)
      NLST2=NLIST(ISYM1,ISYM2,14)
      IF(NLST2.EQ.0) GOTO 100
      VAL2(1)=1.0D00
      VAL2(2)=SQR2
      IX=1
      INCX1=1
      INCX2=NTUV(ISYM1)
      IDIT=1+IOFDIT(ISYM12)
      INCF1=NORB(ISYM12)
      INCF2=1
      IY=1
      INCY1=1
      INCY2=NTGEU(ISYM2)
      CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &            X2(IX),DPT2(IDIT),Y(IY))
      GOTO 100
C  -----------------------------------------------
   2  CONTINUE
C A&BM One-el
      LLST1=LLIST(ISYM1,ISYM2,13)
      NLST1=NLIST(ISYM1,ISYM2,13)
      IF(NLST1.EQ.0) GOTO 92
      VAL1(1)= 3.0D00
      VAL1(2)=-3.0D00
      LLST2=LLIST(ISYM1,ISYM2,15)
      NLST2=NLIST(ISYM1,ISYM2,15)
      IF(NLST2.EQ.0) GOTO 92
      VAL2(1)=-1.0D00
      VAL2(2)= 1.0D00
      IXTI=1
      INCX1=1
      INCX2=NASH(ISYM1)
      IDIT=1+IOFDIT(ISYM12)
      INCF1=NORB(ISYM12)
      INCF2=1
      IY=1
      INCY1=1
      INCY2=NTGTU(ISYM2)
      CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &            X1(IXTI),DPT2(IDIT),Y(IY))
  92  CONTINUE

C A&BM Two-el
      LLST1=LLIST(ISYM1,ISYM2, 4)
      NLST1=NLIST(ISYM1,ISYM2, 4)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)=-1.0D00
      VAL1(2)= 1.0D00
      LLST2=LLIST(ISYM1,ISYM2,15)
      NLST2=NLIST(ISYM1,ISYM2,15)
      IF(NLST2.EQ.0) GOTO 100
      VAL2(1)=-1.0D00
      VAL2(2)= 1.0D00
      IX=1
      INCX1=1
      INCX2=NTUV(ISYM1)
      IDIT=1+IOFDIT(ISYM12)
      INCF1=NORB(ISYM12)
      INCF2=1
      IY=1
      INCY1=1
      INCY2=NTGTU(ISYM2)
      CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &            X2(IX),DPT2(IDIT),Y(IY))
      GOTO 100
C  -----------------------------------------------
   3  CONTINUE

C  A&D  Two-el
      LLST1=LLIST(ISYM1,ISYM2, 1)
      NLST1=NLIST(ISYM1,ISYM2, 1)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)= 1.0D00
      VAL1(2)= 1.0D00
      IX=1
      INCX1=1
      INCX2=NAS1
      IDTA=1+IOFDTA(ISYM12)
      INCF1=1
      INCF2=NORB(ISYM12)
      IY=1+NAS2*IOFCD(ISYM2,ISYM12)
      INCY1=1
      INCY2=NAS2
      INCY3=NAS2*NISH(ISYM1)
      LEN1=NISH(ISYM1)
      LEN2=NSSH(ISYM12)
      CALL MLTMV(IMLTOP,LIST(LLST1),X2(IX),DPT2(IDTA),Y(IY))
      GOTO 100
C  -----------------------------------------------
   4  CONTINUE

C  A&EP One-el
      IF(ISYM2.EQ.ISYM1) THEN
       NT=NASH(ISYM1)
       IF(NT.EQ.0) GOTO 100
       DO ISYMIJ=1,NSYM
        ISYMA=MUL(ISYMIJ,ISYM1)
        NA=NSSH(ISYMA)
        IF(NA.EQ.0) GOTO 94
        LLST1=LLIST(ISYM1,ISYMIJ,14)
        NLST1=NLIST(ISYM1,ISYMIJ,14)
        IF(NLST1.EQ.0) GOTO 94
        VAL1(1)= 1.0D00
        VAL1(2)= SQR2
        IXTI=1
        INCX1=NT
        INCX2=1
        IDIA=1+IOFDIA(ISYMA)
        INCF1=1
        INCF2=NORB(ISYMA)
        IY=1+NAS2*IOFCEP(ISYM2,ISYMA)
        INCY1=NT*NA
        INCY2=1
        INCY3=NT
        LEN1=NT
        LEN2=NA
        CALL MLTMV(IMLTOP,LIST(LLST1),X1(IXTI),DPT2(IDIA),Y(IY))
  94    CONTINUE
       END DO
      END IF
      GOTO 100
C  -----------------------------------------------
   5  CONTINUE

C  A&EM One-el
      IF(ISYM2.EQ.ISYM1) THEN
       NT=NASH(ISYM1)
       IF(NT.EQ.0) GOTO 100
       DO ISYMIJ=1,NSYM
        ISYMA=MUL(ISYMIJ,ISYM1)
        NA=NSSH(ISYMA)
        IF(NA.EQ.0) GOTO 95
        LLST1=LLIST(ISYM1,ISYMIJ,15)
        NLST1=NLIST(ISYM1,ISYMIJ,15)
        IF(NLST1.EQ.0) GOTO 95
        VAL1(1)=-SQR3
        VAL1(2)= SQR3
        IXTI=1
        INCX1=NT
        INCX2=1
        IDIA=1+IOFDIA(ISYMA)
        INCF1=1
        INCF2=NORB(ISYMA)
        IY=1+NAS2*IOFCEM(ISYM2,ISYMA)
        INCY1=NT*NA
        INCY2=1
        INCY3=NT
        LEN1=NT
        LEN2=NA
        CALL MLTMV(IMLTOP,LIST(LLST1),X1(IXTI),DPT2(IDIA),Y(IY))
  95    CONTINUE
       END DO
      END IF
      GOTO 100
C  -----------------------------------------------
   6  CONTINUE

C  BP&EP Two-el
      NA=NSSH(ISYM12)
      IF(NA.EQ.0) GOTO 100
      LLST1=LLIST(ISYM1,ISYM2, 9)
      NLST1=NLIST(ISYM1,ISYM2, 9)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)= SQRI2
      VAL1(2)= SQRI2
      IX=1
      INCX1=1
      INCX2=NAS1
      IDTA=1+IOFDTA(ISYM12)
      INCF1=1
      INCF2=NORB(ISYM12)
      IY=1+NAS2*IOFCEP(ISYM2,ISYM12)
      INCY1=1
      INCY2=NAS2*NA
      INCY3=NAS2
      LEN1=NIS1
      LEN2=NA
      CALL MLTMV(IMLTOP,LIST(LLST1),X2(IX),DPT2(IDTA),Y(IY))
      GOTO 100
C  -----------------------------------------------
   7  CONTINUE

C  BM&EM Two-el
      NA=NSSH(ISYM12)
      IF(NA.EQ.0) GOTO 100
      LLST1=LLIST(ISYM1,ISYM2,10)
      NLST1=NLIST(ISYM1,ISYM2,10)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)=-SQRI6
      VAL1(2)= SQRI6
      IX=1
      INCX1=1
      INCX2=NAS1
      IDTA=1+IOFDTA(ISYM12)
      INCF1=1
      INCF2=NORB(ISYM12)
      IY=1+NAS2*IOFCEM(ISYM2,ISYM12)
      INCY1=1
      INCY2=NAS2*NA
      INCY3=NAS2
      LEN1=NIS1
      LEN2=NA
      CALL MLTMV(IMLTOP,LIST(LLST1),X2(IX),DPT2(IDTA),Y(IY))
      GOTO 100
C  -----------------------------------------------
   8  CONTINUE

C  C&D  One-el
      NI=NISH(ISYM12)
      IF(NI.EQ.0) GOTO 98
      LLST1=LLIST(ISYM1,ISYM2,11)
      NLST1=NLIST(ISYM1,ISYM2,11)
      IF(NLST1.EQ.0) GOTO 98
      VAL1(1)= 2.0D00
      VAL1(2)= 1.0D00
      IXTA=1
      INCX1=1
      INCX2=NASH(ISYM1)
      IDIT=1+IOFDIT(ISYM12)
      INCF1=NORB(ISYM12)
      INCF2=1
      IY=1+NAS2*IOFCD(ISYM2,ISYM1)
      INCY1=1
      INCY2=NAS2*NI
      INCY3=NAS2
      LEN1=NSSH(ISYM1)
      LEN2=NI
      CALL MLTMV(IMLTOP,LIST(LLST1),X1(IXTA),DPT2(IDIT),Y(IY))
  98  CONTINUE

C  C&D  Two-el
      NI=NISH(ISYM12)
      IF(NI.EQ.0) GOTO 100
      LLST1=LLIST(ISYM1,ISYM2, 2)
      NLST1=NLIST(ISYM1,ISYM2, 2)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)=-1.0D00
      VAL1(2)=-1.0D00
      IX=1
      INCX1=1
      INCX2=NAS1
      IDIT=1+IOFDIT(ISYM12)
      INCF1=NORB(ISYM12)
      INCF2=1
      IY=1+NAS2*IOFCD(ISYM2,ISYM1)
      INCY1=1
      INCY2=NAS2*NI
      INCY3=NAS2
      LEN1=NSSH(ISYM1)
      LEN2=NI
      CALL MLTMV(IMLTOP,LIST(LLST1),X2(IX),DPT2(IDIT),Y(IY))
      GOTO 100
C  -----------------------------------------------
   9  CONTINUE

C  C&FP One-el
      LLST1=LLIST(ISYM1,ISYM2,12)
      NLST1=NLIST(ISYM1,ISYM2,12)
      IF(NLST1.EQ.0) GOTO 99
      VAL1(1)=-1.0D00
      VAL1(2)=-2.0D00
      LLST2=LLIST(ISYM1,ISYM2,16)
      NLST2=NLIST(ISYM1,ISYM2,16)
      IF(NLST2.EQ.0) GOTO 99
      VAL2(1)=1.0D00
      VAL2(2)=SQR2
      IXTA=1
      INCX1=1
      INCX2=NASH(ISYM1)
      IDTA=1+IOFDTA(ISYM12)
      INCF1=1
      INCF2=NORB(ISYM12)
      IY=1
      INCY1=1
      INCY2=NTGEU(ISYM2)
      CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &            X1(IXTA),DPT2(IDTA),Y(IY))
  99  CONTINUE

C  C&FP Two-el
      LLST1=LLIST(ISYM1,ISYM2, 5)
      NLST1=NLIST(ISYM1,ISYM2, 5)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)= 1.0D00
      VAL1(2)= 2.0D00
      LLST2=LLIST(ISYM1,ISYM2,16)
      NLST2=NLIST(ISYM1,ISYM2,16)
      IF(NLST2.EQ.0) GOTO 100
      VAL2(1)=1.0D00
      VAL2(2)=SQR2
      IX=1
      INCX1=1
      INCX2=NTUV(ISYM1)
      IDTA=1+IOFDTA(ISYM12)
      INCF1=1
      INCF2=NORB(ISYM12)
      IY=1
      INCY1=1
      INCY2=NTGEU(ISYM2)
      CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &            X2(IX),DPT2(IDTA),Y(IY))
      GOTO 100
C  -----------------------------------------------
  10  CONTINUE

C  C&FM One-el
      LLST1=LLIST(ISYM1,ISYM2,13)
      NLST1=NLIST(ISYM1,ISYM2,13)
      IF(NLST1.EQ.0) GOTO 910
      VAL1(1)=-1.0D00
      VAL1(2)= 1.0D00
      LLST2=LLIST(ISYM1,ISYM2,17)
      NLST2=NLIST(ISYM1,ISYM2,17)
      IF(NLST2.EQ.0) GOTO 910
      VAL2(1)=1.0D00
      VAL2(2)=-1.0D00
      IXTA=1
      INCX1=1
      INCX2=NASH(ISYM1)
      IDTA=1+IOFDTA(ISYM12)
      INCF1=1
      INCF2=NORB(ISYM12)
      IY=1
      INCY1=1
      INCY2=NTGTU(ISYM2)
      CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &            X1(IXTA),DPT2(IDTA),Y(IY))
 910  CONTINUE

C  C&FM Two-el
      LLST1=LLIST(ISYM1,ISYM2, 6)
      NLST1=NLIST(ISYM1,ISYM2, 6)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)=-1.0D00
      VAL1(2)= 1.0D00
      LLST2=LLIST(ISYM1,ISYM2,17)
      NLST2=NLIST(ISYM1,ISYM2,17)
      IF(NLST2.EQ.0) GOTO 100
      VAL2(1)=1.0D00
      VAL2(2)=-1.0D00
      IX=1
      INCX1=1
      INCX2=NTUV(ISYM1)
      IDTA=1+IOFDTA(ISYM12)
      INCF1=1
      INCF2=NORB(ISYM12)
      IY=1
      INCY1=1
      INCY2=NTGTU(ISYM2)
      CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &            X2(IX),DPT2(IDTA),Y(IY))
      GOTO 100
C  -----------------------------------------------
  11  CONTINUE

C  C&GP One-el
      IF(ISYM2.EQ.ISYM1) THEN
       NT=NASH(ISYM1)
       IF(NT.EQ.0) GOTO 100
       DO ISYMAB=1,NSYM
        ISYMI=MUL(ISYMAB,ISYM1)
        NI=NISH(ISYMI)
        IF(NI.EQ.0) GOTO 911
        LLST1=LLIST(ISYM1,ISYMAB,16)
        NLST1=NLIST(ISYM1,ISYMAB,16)
        IF(NLST1.EQ.0) GOTO 911
        VAL1(1)= SQRI2
        VAL1(2)= 1.0D00
        IXTA=1
        INCX1=NT
        INCX2=1
        IDIA=1+IOFDIA(ISYMI)
        INCF1=NORB(ISYMI)
        INCF2=1
        IY=1+NAS2*IOFCGP(ISYM2,ISYMI)
        INCY1=NT*NI
        INCY2=1
        INCY3=NT
        LEN1=NT
        LEN2=NI
        CALL MLTMV(IMLTOP,LIST(LLST1),X1(IXTA),DPT2(IDIA),Y(IY))
 911    CONTINUE
       END DO
      END IF
      GOTO 100
C  -----------------------------------------------
  12  CONTINUE

C  C&GM One-el
      IF(ISYM2.EQ.ISYM1) THEN
       NT=NASH(ISYM1)
       IF(NT.EQ.0) GOTO 100
       DO ISYMAB=1,NSYM
        ISYMI=MUL(ISYMAB,ISYM1)
        NI=NISH(ISYMI)
        IF(NI.EQ.0) GOTO 912
        LLST1=LLIST(ISYM1,ISYMAB,17)
        NLST1=NLIST(ISYM1,ISYMAB,17)
        IF(NLST1.EQ.0) GOTO 912
        VAL1(1)= SQR32
        VAL1(2)=-SQR32
        IXTA=1
        INCX1=NT
        INCX2=1
        IDIA=1+IOFDIA(ISYMI)
        INCF1=NORB(ISYMI)
        INCF2=1
        IY=1+NAS2*IOFCGM(ISYM2,ISYMI)
        INCY1=NT*NI
        INCY2=1
        INCY3=NT
        LEN1=NT
        LEN2=NI
        CALL MLTMV(IMLTOP,LIST(LLST1),X1(IXTA),DPT2(IDIA),Y(IY))
 912    CONTINUE
       END DO
      END IF
      GOTO 100
C  -----------------------------------------------
  13  CONTINUE

C  D&EP One-el
      IF(ISYM1.EQ.1) THEN
       NT=NASH(ISYM2)
       IF(NT.EQ.0) GOTO 913
CPAM98 Added line:
       IOXIA=0
       DO ISYMI=1,NSYM
        NI=NISH(ISYMI)
        IF(NI.EQ.0) GOTO 813
        NA=NSSH(ISYMI)
        IF(NA.EQ.0) GOTO 813
        ISYMIJ=MUL(ISYMI,ISYM2)
        LLST1=LLIST(ISYMI,ISYMIJ,14)
        NLST1=NLIST(ISYMI,ISYMIJ,14)
        IF(NLST1.EQ.0) GOTO 813
        VAL1(1)= SQRI2
        VAL1(2)= 1.0D00
CPAM98        IXIA=1+IOFDIA(ISYMI)
CPAM98 Added line:
        IXIA=IOXIA+1
        INCX1=1
        INCX2=NI
        IDIT=1+IOFDIT(ISYM2)
        INCF1=1
        INCF2=NORB(ISYM2)
        IY=1+NAS2*IOFCEP(ISYM2,ISYMI)
        INCY1=NT*NA
        INCY2=NT
        INCY3=1
        LEN1=NA
        LEN2=NT
        CALL MLTMV(IMLTOP,LIST(LLST1),X1(IXIA),DPT2(IDIT),Y(IY))
 813    CONTINUE
CPAM98 Added line:
        IOXIA=IOXIA+NI*NA
       END DO
 913   CONTINUE
      END IF

C  D&EP Two-el
      LLST1=LLIST(ISYM1,ISYM2, 7)
      NLST1=NLIST(ISYM1,ISYM2, 7)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)=-1.0D00
      VAL1(2)=-1.0D00
      NU=NASH(ISYM2)
      IF(NU.EQ.0) GOTO 100
      DO ISYMA=1,NSYM
        NA=NSSH(ISYMA)
        IF(NA.EQ.0) GOTO 713
        ISYMI=MUL(ISYMA,ISYM1)
        NI=NISH(ISYMI)
        IF(NI.EQ.0) GOTO 713
        ISYMIJ=MUL(ISYMI,ISYM12)
        LLST2=LLIST(ISYMI,ISYMIJ,14)
        NLST2=NLIST(ISYMI,ISYMIJ,14)
        IF(NLST2.EQ.0) GOTO 713
        VAL2(1)= SQRI2
        VAL2(2)= 1.0D00
        IX=1+NAS1*IOFCD(ISYM1,ISYMA)
        INCX1=1
        INCX2=NAS1
        INCX3=NAS1*NI
        IDIT=1+IOFDIT(ISYM12)
        INCF1=NORB(ISYM12)
        INCF2=1
        IY=1+NU*IOFCEP(ISYM2,ISYMA)
        INCY1=1
        INCY2=NU*NA
        INCY3=NU
        LEN1=NA
        CALL MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),
     &              X2(IX),DPT2(IDIT),Y(IY))
 713    CONTINUE
      END DO
      GOTO 100
C  -----------------------------------------------
  14  CONTINUE

C  D&EM One-el
      IF(ISYM1.EQ.1) THEN
       NT=NASH(ISYM2)
       IF(NT.EQ.0) GOTO 914
CPAM98 Added line:
       IOXIA=0
       DO ISYMI=1,NSYM
        NI=NISH(ISYMI)
        IF(NI.EQ.0) GOTO 814
        NA=NSSH(ISYMI)
        IF(NI.EQ.0) GOTO 814
        ISYMIJ=MUL(ISYMI,ISYM2)
        LLST1=LLIST(ISYMI,ISYMIJ,15)
        NLST1=NLIST(ISYMI,ISYMIJ,15)
        IF(NLST1.EQ.0) GOTO 814
        VAL1(1)= SQR32
        VAL1(2)=-SQR32
CPAM98        IXIA=1+IOFDIA(ISYMI)
CPAM98 Added line:
        IXIA=IOXIA+1
        INCX1=1
        INCX2=NI
        IDIT=1+IOFDIT(ISYM2)
        INCF1=1
        INCF2=NORB(ISYM2)
        IY=1+NAS2*IOFCEM(ISYM2,ISYMI)
        INCY1=NT*NA
        INCY2=NT
        INCY3=1
        LEN1=NA
        LEN2=NT
        CALL MLTMV(IMLTOP,LIST(LLST1),X1(IXIA),DPT2(IDIT),Y(IY))
 814    CONTINUE
CPAM98 Added line:
        IOXIA=IOXIA+NI*NA
       END DO
 914   CONTINUE
      END IF

C  D&EM Two-el
      LLST1=LLIST(ISYM1,ISYM2, 7)
      NLST1=NLIST(ISYM1,ISYM2, 7)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)=-1.0D00
      VAL1(2)= 1.0D00
      NU=NASH(ISYM2)
      IF(NU.EQ.0) GOTO 100
      DO ISYMA=1,NSYM
        NA=NSSH(ISYMA)
        IF(NA.EQ.0) GOTO 714
        ISYMI=MUL(ISYMA,ISYM1)
        NI=NISH(ISYMI)
        IF(NI.EQ.0) GOTO 714
        ISYMIJ=MUL(ISYMI,ISYM12)
        LLST2=LLIST(ISYMI,ISYMIJ,15)
        NLST2=NLIST(ISYMI,ISYMIJ,15)
        IF(NLST2.EQ.0) GOTO 714
        VAL2(1)= SQRI6
        VAL2(2)=-SQRI6
        IX=1+NAS1*IOFCD(ISYM1,ISYMA)
        INCX1=1
        INCX2=NAS1
        INCX3=NAS1*NI
        IDIT=1+IOFDIT(ISYM12)
        INCF1=NORB(ISYM12)
        INCF2=1
        IY=1+NU*IOFCEM(ISYM2,ISYMA)
        INCY1=1
        INCY2=NU*NA
        INCY3=NU
        LEN1=NA
        CALL MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),
     &              X2(IX),DPT2(IDIT),Y(IY))
 714    CONTINUE
      END DO
      GOTO 100
C  -----------------------------------------------
  15  CONTINUE

C  D&GP Two-el
      LLST1=LLIST(ISYM1,ISYM2, 8)
      NLST1=NLIST(ISYM1,ISYM2, 8)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)= 1.0D00
      VAL1(2)= 1.0D00
      NU=NASH(ISYM2)
      IF(NU.EQ.0) GOTO 100
      DO ISYMA=1,NSYM
        ISYMI=MUL(ISYMA,ISYM1)
        NI=NISH(ISYMI)
        IF(NI.EQ.0) GOTO 915
        ISYMAB=MUL(ISYMA,ISYM12)
        LLST2=LLIST(ISYMA,ISYMAB,16)
        NLST2=NLIST(ISYMA,ISYMAB,16)
        IF(NLST2.EQ.0) GOTO 915
        VAL2(1)= SQRI2
        VAL2(2)= 1.0D00
        IX=1+NAS1*IOFCD(ISYM1,ISYMA)
        INCX1=1
        INCX2=NAS1*NI
        INCX3=NAS1
        IDTA=1+IOFDTA(ISYM12)
        INCF1=1
        INCF2=NORB(ISYM12)
        IY=1+NU*IOFCGP(ISYM2,ISYMI)
        INCY1=1
        INCY2=NU*NI
        INCY3=NU
        LEN1=NI
        CALL MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),
     &              X2(IX),DPT2(IDTA),Y(IY))
 915    CONTINUE
      END DO
      GOTO 100
C  -----------------------------------------------
  16  CONTINUE

C  D&GM Two-el
      LLST1=LLIST(ISYM1,ISYM2, 8)
      NLST1=NLIST(ISYM1,ISYM2, 8)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)=-1.0D00
      VAL1(2)= 1.0D00
      NU=NASH(ISYM2)
      IF(NU.EQ.0) GOTO 100
      DO ISYMA=1,NSYM
        ISYMI=MUL(ISYMA,ISYM1)
        NI=NISH(ISYMI)
        IF(NI.EQ.0) GOTO 916
        ISYMAB=MUL(ISYMA,ISYM12)
        LLST2=LLIST(ISYMA,ISYMAB,17)
        NLST2=NLIST(ISYMA,ISYMAB,17)
        IF(NLST2.EQ.0) GOTO 916
        VAL2(1)= SQRI6
        VAL2(2)=-SQRI6
        IX=1+NAS1*IOFCD(ISYM1,ISYMA)
        INCX1=1
        INCX2=NAS1*NI
        INCX3=NAS1
        IDTA=1+IOFDTA(ISYM12)
        INCF1=1
        INCF2=NORB(ISYM12)
        IY=1+NU*IOFCGM(ISYM2,ISYMI)
        INCY1=1
        INCY2=NU*NI
        INCY3=NU
        LEN1=NI
        CALL MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),
     &              X2(IX),DPT2(IDTA),Y(IY))
 916    CONTINUE
      END DO
      GOTO 100
C  -----------------------------------------------
  17  CONTINUE

C  EP&HP Two-el
      NA=NSSH(ISYM12)
      IF(NA.EQ.0) GOTO 100
      LLST1=LLIST(ISYM12,ISYM2,16)
      NLST1=NLIST(ISYM12,ISYM2,16)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)= SQRI2
      VAL1(2)= 1.0D00
      IX=1+NAS1*IOFCEP(ISYM1,ISYM12)
      INCX1=NAS1
      INCX2=1
      INCX3=NAS1*NA
      IDTA=1+IOFDTA(ISYM1)
      INCF1=NORB(ISYM1)
      INCF2=1
      IY=1
      INCY1=1
      INCY2=NAS2
      LEN1=NAS1
      LEN2=NIGEJ(ISYM2)
      CALL MLTR1(IMLTOP,LIST(LLST1),X2(IX),DPT2(IDTA),Y(IY))
      GOTO 100
C  -----------------------------------------------
  18  CONTINUE

C  EM&HM Two-el
      NA=NSSH(ISYM12)
      IF(NA.EQ.0) GOTO 100
      LLST1=LLIST(ISYM12,ISYM2,17)
      NLST1=NLIST(ISYM12,ISYM2,17)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)= SQRI2
      VAL1(2)=-SQRI2
      IX=1+NAS1*IOFCEM(ISYM1,ISYM12)
      INCX1=NAS1
      INCX2=1
      INCX3=NAS1*NA
      IDTA=1+IOFDTA(ISYM1)
      INCF1=NORB(ISYM1)
      INCF2=1
      IY=1
      INCY1=1
      INCY2=NAS2
      LEN1=NAS1
      LEN2=NIGTJ(ISYM2)
      CALL MLTR1(IMLTOP,LIST(LLST1),X2(IX),DPT2(IDTA),Y(IY))
      GOTO 100
C  -----------------------------------------------
  19  CONTINUE

C  FP&GP Two-el
      NI=NISH(ISYM12)
      IF(NI.EQ.0) GOTO 100
      LLST1=LLIST(ISYM1,ISYM2, 9)
      NLST1=NLIST(ISYM1,ISYM2, 9)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)=-SQRI2
      VAL1(2)=-SQRI2
      IX=1
      INCX1=1
      INCX2=NAS1
      IDIT=1+IOFDIT(ISYM12)
      INCF1=NORB(ISYM12)
      INCF2=1
      IY=1+NAS2*IOFCGP(ISYM2,ISYM12)
      INCY1=1
      INCY2=NAS2*NI
      INCY3=NAS2
      LEN1=NIS1
      LEN2=NI
      CALL MLTMV(IMLTOP,LIST(LLST1),X2(IX),DPT2(IDIT),Y(IY))
      GOTO 100
C  -----------------------------------------------
  20  CONTINUE

C  FM&GM Two-el
      NI=NISH(ISYM12)
      IF(NI.EQ.0) GOTO 100
      LLST1=LLIST(ISYM1,ISYM2,10)
      NLST1=NLIST(ISYM1,ISYM2,10)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)=-SQRI6
      VAL1(2)= SQRI6
      IX=1
      INCX1=1
      INCX2=NAS1
      IDIT=1+IOFDIT(ISYM12)
      INCF1=NORB(ISYM12)
      INCF2=1
      IY=1+NAS2*IOFCGM(ISYM2,ISYM12)
      INCY1=1
      INCY2=NAS2*NI
      INCY3=NAS2
      LEN1=NIS1
      LEN2=NI
      CALL MLTMV(IMLTOP,LIST(LLST1),X2(IX),DPT2(IDIT),Y(IY))
      GOTO 100
C  -----------------------------------------------
  21  CONTINUE

C  GP&HP Two-el
      LLST1=LLIST(ISYM12,ISYM2,14)
      NLST1=NLIST(ISYM12,ISYM2,14)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)= -SQRI2
      VAL1(2)= -1.0D00
      IX=1+NAS1*IOFCGP(ISYM1,ISYM12)
      INCX1=NAS1
      INCX2=1
      INCX3=NAS1*NISH(ISYM12)
      IDIT=1+IOFDIT(ISYM1)
      INCF1=1
      INCF2=NORB(ISYM1)
      IY=1
      INCY1=NAS2
      INCY2=1
      LEN1=NAS1
      LEN2=NAS2
      CALL MLTR1(IMLTOP,LIST(LLST1),X2(IX),DPT2(IDIT),Y(IY))
      GOTO 100
C  -----------------------------------------------
  22  CONTINUE

C  GM&HM Two-el
      LLST1=LLIST(ISYM12,ISYM2,15)
      NLST1=NLIST(ISYM12,ISYM2,15)
      IF(NLST1.EQ.0) GOTO 100
      VAL1(1)= SQRI2
      VAL1(2)=-SQRI2
      IX=1+NAS1*IOFCGM(ISYM1,ISYM12)
      INCX1=NAS1
      INCX2=1
      INCX3=NAS1*NISH(ISYM12)
      IDIT=1+IOFDIT(ISYM1)
      INCF1=1
      INCF2=NORB(ISYM1)
      IY=1
      INCY1=NAS2
      INCY2=1
      LEN1=NAS1
      LEN2=NAS2
      CALL MLTR1(IMLTOP,LIST(LLST1),X2(IX),DPT2(IDIT),Y(IY))
      GOTO 100
C  -----------------------------------------------
  23  CONTINUE

C  D&HP One-el
      IF(ISYM1.NE.1) GOTO 100
      IOXIA=0
      DO ISYMI=1,NSYM
       NI=NISH(ISYMI)
       IF(NI.EQ.0) GOTO 923
       ISYMA=ISYMI
       NA=NSSH(ISYMA)
       IF(NA.EQ.0) GOTO 923
       ISYMJ=MUL(ISYMI,ISYM2)
       NJ=NISH(ISYMJ)
       IF(NJ.EQ.0) GOTO 923
       ISYMB=ISYMJ
       NB=NSSH(ISYMB)
       IF(NB.EQ.0) GOTO 923
       LLST1=LLIST(ISYMI,ISYM2,14)
       NLST1=NLIST(ISYMI,ISYM2,14)
       IF(NLST1.EQ.0) GOTO 923
       VAL1(1)=SQRI2
       VAL1(2)=1.0D0
       LLST2=LLIST(ISYMA,ISYM2,16)
       NLST2=NLIST(ISYMA,ISYM2,16)
       IF(NLST2.EQ.0) GOTO 923
       VAL2(1)=SQRI2
       VAL2(2)=1.0D0
       IXIA=IOXIA+1
       INCX1=1
       INCX2=NI
       IDJB=1+IOFDIA(ISYMB)
       INCF1=1
       INCF2=NORB(ISYMJ)
       IY=1
       INCY1=NAGEB(ISYM2)
       INCY2=1
       CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &              X1(IXIA),DPT2(IDJB),Y(IY))
923    CONTINUE
       IOXIA=IOXIA+NI*NA
      END DO
      GOTO 100

C  ---------------------------

  24  CONTINUE
C  D&HM One-el
      IF(ISYM1.NE.1) GOTO 100
      IOXIA=0
      DO ISYMI=1,NSYM
       NI=NISH(ISYMI)
       IF(NI.EQ.0) GOTO 924
       ISYMA=ISYMI
       NA=NSSH(ISYMA)
       IF(NA.EQ.0) GOTO 924
       ISYMJ=MUL(ISYMI,ISYM2)
       NJ=NISH(ISYMJ)
       IF(NJ.EQ.0) GOTO 924
       ISYMB=ISYMJ
       NB=NSSH(ISYMB)
       IF(NB.EQ.0) GOTO 924
       LLST1=LLIST(ISYMI,ISYM2,15)
       NLST1=NLIST(ISYMI,ISYM2,15)
       IF(NLST1.EQ.0) GOTO 924
       VAL1(1)=SQR3*0.5D0
       VAL1(2)=-VAL1(1)
       LLST2=LLIST(ISYMA,ISYM2,17)
       NLST2=NLIST(ISYMA,ISYM2,17)
       IF(NLST2.EQ.0) GOTO 924
       VAL2(1)=1.0D0
       VAL2(2)=-1.0D0
       IXIA=IOXIA+1
       INCX1=1
       INCX2=NI
       IDJB=1+IOFDIA(ISYMB)
       INCF1=1
       INCF2=NORB(ISYMJ)
       IY=1
       INCY1=NAGTB(ISYM2)
       INCY2=1
       CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &            X1(IXIA),DPT2(IDJB),Y(IY))
 924   CONTINUE
       IOXIA=IOXIA+NI*NA
      END DO
      GOTO 100
C  -----------------------------------------------
 100  CONTINUE

      RETURN
      END
