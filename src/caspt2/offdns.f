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
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use constants, only: Half, One, Two, Three, Six
      use EQSOLV, only: IFCOUP, NLIST, LLIST
      use Sigma_data, only: VAL1, VAL2, IFTEST, INCF1, INCF2, INCX1,
     &                      INCX2, INCX3, INCY1, INCY2, INCY3, LEN1,
     &                      LEN2, NLST1, NLST2
      use caspt2_module, only: NSYM, NISUP, NASUP, NASH, NORB,
     &                         NSSH, NISH, NIGEJ, NIGTJ, NAGEB, NAGTB,
     &                         NTGEU, NTUV, NTGTU
      IMPLICIT None

      integer(kind=iwp), intent(in)::  ISYM1,ICASE1,ISYM2,ICASE2
      real(kind=wp), intent(inout):: X1(*),X2(*),Y(*)
      real(kind=wp), intent(inout):: DPT2(*)
      integer(kind=iwp), intent(in):: LIST(*)

      integer(kind=iwp) IOFDIT(8),IOFDIA(8),IOFDTA(8)
      integer(kind=iwp) IOFCD(8,8),IOFCEP(8,8),IOFCEM(8,8),IOFCGP(8,8),
     &                  IOFCGM(8,8)
C Various constants:
      real(kind=wp), parameter:: SQR2=SQRT(Two), SQR3=SQRT(Three),
     &                           SQR6=SQRT(Six), SQRI2=One/SQR2,
     &                           SQRI6=One/SQR6, SQR32=SQR3*SQRI2
      Integer(kind=iwp) ICD, ICEM, ICEP, ICGM, ICGP, IDIA, IDIT, IDJB,
     &                  IDOFF, IDTA, IJSYM, IMLTOP, IOXIA, ISYM, ISYM12,
     &                  ISYMA, ISYMAB, ISYMB, ISYMI, ISYMIJ, ISYMJ, IX,
     &                  IXIA, IXTA, IXTI, IY, JSYM, KOD, LLST1, LLST2,
     &                  NA, NAS1, NAS2, NB, NI, NIS1, NJ, NO, NT, NU
C Compute off-diagonal contributions to a trans density matrix.
C Sub-diagonal blocks only. If a density matrix is required,
C i.e. both wave functions equal, use symmetry. Else, call
C again with wave functions interchanged.


      NA = 0 ! dummy initialize
*
      IFTEST=0
      IMLTOP=2
      KOD=IFCOUP(ICASE2,ICASE1)
      IF (KOD==0) RETURN

      ISYM12=Mul(ISYM1,ISYM2)
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
          IJSYM=Mul(ISYM,JSYM)
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

      Select case (KOD)
C  -----------------------------------------------
      Case(1)
C  A&BP One-el
      NLST1=NLIST(ISYM1,ISYM2,12)
      NLST2=NLIST(ISYM1,ISYM2,14)
      IF(NLST1*NLST2/=0) THEN
        LLST1=LLIST(ISYM1,ISYM2,12)
        VAL1(1)=One
        VAL1(2)=Two
        LLST2=LLIST(ISYM1,ISYM2,14)
        VAL2(1)=One
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
     &              X1(IXTI),DPT2(IDIT),Y(IY))
      END IF

C  A&BP Two-el
      NLST1=NLIST(ISYM1,ISYM2, 3)
      NLST2=NLIST(ISYM1,ISYM2,14)
      IF(NLST1*NLST2/=0) THEN
        LLST1=LLIST(ISYM1,ISYM2, 3)
        VAL1(1)=-One
        VAL1(2)=-Two
        LLST2=LLIST(ISYM1,ISYM2,14)
        VAL2(1)=One
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
     &              X2(IX),DPT2(IDIT),Y(IY))
      END IF
C  -----------------------------------------------
      Case(2)
C A&BM One-el
      NLST1=NLIST(ISYM1,ISYM2,13)
      NLST2=NLIST(ISYM1,ISYM2,15)
      IF(NLST1*NLST2/=0) THEN
        LLST1=LLIST(ISYM1,ISYM2,13)
        VAL1(1)= Three
        VAL1(2)=-Three
        LLST2=LLIST(ISYM1,ISYM2,15)
* Original:
*       VAL2(1)=-One
*       VAL2(2)= One
* Fix for sign error noted by Takeshi, May 2015:
        VAL2(1)= One
        VAL2(2)=-One
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
     &              X1(IXTI),DPT2(IDIT),Y(IY))
      END IF

C A&BM Two-el
      NLST1=NLIST(ISYM1,ISYM2, 4)
      NLST2=NLIST(ISYM1,ISYM2,15)
      IF(NLST1*NLST2/=0) THEN
        LLST1=LLIST(ISYM1,ISYM2, 4)
        VAL1(1)=-One
        VAL1(2)= One
        LLST2=LLIST(ISYM1,ISYM2,15)
* Original:
*       VAL2(1)=-One
*       VAL2(2)= One
* Fix for sign error noted by Takeshi, May 2015:
        VAL2(1)= One
        VAL2(2)=-One
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
     &              X2(IX),DPT2(IDIT),Y(IY))
      END IF
C  -----------------------------------------------
      Case(3)

C  A&D  Two-el
      LLST1=LLIST(ISYM1,ISYM2, 1)
      NLST1=NLIST(ISYM1,ISYM2, 1)
      IF(NLST1/=0) THEN
        VAL1(1)= One
        VAL1(2)= One
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
      END IF
C  -----------------------------------------------
      Case(4)

C  A&EP One-el
      NT=NASH(ISYM1)
      IF(ISYM2==ISYM1 .AND. NT/=0) THEN
       DO ISYMIJ=1,NSYM
        ISYMA=Mul(ISYMIJ,ISYM1)
        NA=NSSH(ISYMA)
        NLST1=NLIST(ISYM1,ISYMIJ,14)
        IF(NA*NLST1/=0) THEN
          LLST1=LLIST(ISYM1,ISYMIJ,14)
          VAL1(1)= One
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
        END IF
       END DO
      END IF
C  -----------------------------------------------
      Case(5)

C  A&EM One-el
      NT=NASH(ISYM1)
      IF(ISYM2==ISYM1 .AND. NT/=0) THEN
       DO ISYMIJ=1,NSYM
        ISYMA=Mul(ISYMIJ,ISYM1)
        NA=NSSH(ISYMA)
        NLST1=NLIST(ISYM1,ISYMIJ,15)
        IF(NA*NLST1/=0) THEN
          LLST1=LLIST(ISYM1,ISYMIJ,15)
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
        END IF
       END DO
      END IF
C  -----------------------------------------------
      Case(6)

C  BP&EP Two-el
      NA=NSSH(ISYM12)
      NLST1=NLIST(ISYM1,ISYM2, 9)
      IF(NA*NLST1/=0) THEN
        LLST1=LLIST(ISYM1,ISYM2, 9)
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
      END IF
C  -----------------------------------------------
      Case(7)

C  BM&EM Two-el
      NA=NSSH(ISYM12)
      NLST1=NLIST(ISYM1,ISYM2,10)
      IF(NA*NLST1/=0) THEN
        LLST1=LLIST(ISYM1,ISYM2,10)
* Original:
*       VAL1(1)=-SQRI6
*       VAL1(2)= SQRI6
* Fix for sign error noted by Takeshi, May 2015:
        VAL1(1)= SQRI6
        VAL1(2)=-SQRI6
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
      END IF
C  -----------------------------------------------
      Case(8)

C  C&D  One-el
      NI=NISH(ISYM12)
      NLST1=NLIST(ISYM1,ISYM2,11)
      IF(NI*NLST1/=0) THEN
        LLST1=LLIST(ISYM1,ISYM2,11)
        VAL1(1)= Two
        VAL1(2)= One
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
      END IF

C  C&D  Two-el
      NI=NISH(ISYM12)
      NLST1=NLIST(ISYM1,ISYM2, 2)
      IF(NI*NLST1/=0) THEN
        LLST1=LLIST(ISYM1,ISYM2, 2)
        VAL1(1)=-One
        VAL1(2)=-One
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
      END IF
C  -----------------------------------------------
      Case(9)

C  C&FP One-el
      NLST1=NLIST(ISYM1,ISYM2,12)
      NLST2=NLIST(ISYM1,ISYM2,16)
      IF(NLST1/=0) THEN
        LLST1=LLIST(ISYM1,ISYM2,12)
        VAL1(1)=-One
        VAL1(2)=-Two
        LLST2=LLIST(ISYM1,ISYM2,16)
        VAL2(1)=One
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
     &              X1(IXTA),DPT2(IDTA),Y(IY))
      END IF

C  C&FP Two-el
      NLST1=NLIST(ISYM1,ISYM2, 5)
      NLST2=NLIST(ISYM1,ISYM2,16)
      IF(NLST1*NLST2/=0) THEN
        LLST1=LLIST(ISYM1,ISYM2, 5)
        VAL1(1)= One
        VAL1(2)= Two
        LLST2=LLIST(ISYM1,ISYM2,16)
        VAL2(1)=One
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
     &              X2(IX),DPT2(IDTA),Y(IY))
      END IF
C  -----------------------------------------------
      Case(10)

C  C&FM One-el
      NLST1=NLIST(ISYM1,ISYM2,13)
      NLST2=NLIST(ISYM1,ISYM2,17)
      IF(NLST1*NLST2/=0) THEN
        LLST1=LLIST(ISYM1,ISYM2,13)
        VAL1(1)=-One
        VAL1(2)= One
        LLST2=LLIST(ISYM1,ISYM2,17)
        VAL2(1)=One
        VAL2(2)=-One
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
     &              X1(IXTA),DPT2(IDTA),Y(IY))
      END IF

C  C&FM Two-el
      NLST1=NLIST(ISYM1,ISYM2, 6)
      NLST2=NLIST(ISYM1,ISYM2,17)
      IF(NLST1*NLST2/=0) THEN
        LLST1=LLIST(ISYM1,ISYM2, 6)
        VAL1(1)=-One
        VAL1(2)= One
        LLST2=LLIST(ISYM1,ISYM2,17)
        VAL2(1)=One
        VAL2(2)=-One
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
     &              X2(IX),DPT2(IDTA),Y(IY))
      END IF
C  -----------------------------------------------
      Case(11)

C  C&GP One-el
      NT=NASH(ISYM1)
      IF(ISYM2.EQ.ISYM1 .AND. NT/=0) THEN
       DO ISYMAB=1,NSYM
        ISYMI=Mul(ISYMAB,ISYM1)
        NI=NISH(ISYMI)
        NLST1=NLIST(ISYM1,ISYMAB,16)
        IF(NI*NLST1/=0) THEN
          LLST1=LLIST(ISYM1,ISYMAB,16)
          VAL1(1)= SQRI2
          VAL1(2)= One
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
        END IF
       END DO
      END IF
C  -----------------------------------------------
      Case(12)

C  C&GM One-el
      NT=NASH(ISYM1)
      IF(ISYM2.EQ.ISYM1 .AND. NT/=0) THEN
       DO ISYMAB=1,NSYM
        ISYMI=Mul(ISYMAB,ISYM1)
        NI=NISH(ISYMI)
        NLST1=NLIST(ISYM1,ISYMAB,17)
        IF(NI*NLST1/=0) THEN
          LLST1=LLIST(ISYM1,ISYMAB,17)
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
        END IF
       END DO
      END IF
C  -----------------------------------------------
      Case(13)

C  D&EP One-el
      NT=NASH(ISYM2)
      IF(ISYM1==1 .AND. NT/=0) THEN
       IOXIA=0
       DO ISYMI=1,NSYM
        NI=NISH(ISYMI)
        NA=NSSH(ISYMI)
        ISYMIJ=Mul(ISYMI,ISYM2)
        LLST1=LLIST(ISYMI,ISYMIJ,14)
        NLST1=NLIST(ISYMI,ISYMIJ,14)
        IF(NI*NA*NLST1/=0) THEN
          VAL1(1)= SQRI2
          VAL1(2)= One
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
        END IF
        IOXIA=IOXIA+NI*NA
       END DO
      END IF

C  D&EP Two-el
      NLST1=NLIST(ISYM1,ISYM2, 7)
      NU=NASH(ISYM2)
      IF(NLST1*NU/=0) THEN
      LLST1=LLIST(ISYM1,ISYM2, 7)
      VAL1(1)=-One
      VAL1(2)=-One
      DO ISYMA=1,NSYM
        NA=NSSH(ISYMA)
        ISYMI=Mul(ISYMA,ISYM1)
        NI=NISH(ISYMI)
        ISYMIJ=Mul(ISYMI,ISYM12)
        LLST2=LLIST(ISYMI,ISYMIJ,14)
        NLST2=NLIST(ISYMI,ISYMIJ,14)
        IF(NA*NI*NLST2/=0) THEN
          VAL2(1)= SQRI2
          VAL2(2)= One
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
     &                X2(IX),DPT2(IDIT),Y(IY))
        END IF
      END DO
      END IF
C  -----------------------------------------------
      Case(14)

C  D&EM One-el
      NT=NASH(ISYM2)
      IF(ISYM1==1 .AND. NT/=0) THEN
       IOXIA=0
       DO ISYMI=1,NSYM
        NI=NISH(ISYMI)
        NA=NSSH(ISYMI)
        ISYMIJ=Mul(ISYMI,ISYM2)
        LLST1=LLIST(ISYMI,ISYMIJ,15)
        NLST1=NLIST(ISYMI,ISYMIJ,15)
        IF(NI*NLST1/=0) THEN
          VAL1(1)= SQR32
          VAL1(2)=-SQR32
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
        END IF
        IOXIA=IOXIA+NI*NA
       END DO
      END IF

C  D&EM Two-el
      NLST1=NLIST(ISYM1,ISYM2, 7)
      NU=NASH(ISYM2)
      IF(NLST1*NU/=0) THEN
      LLST1=LLIST(ISYM1,ISYM2, 7)
      VAL1(1)=-One
      VAL1(2)= One
      DO ISYMA=1,NSYM
        NA=NSSH(ISYMA)
        ISYMI=Mul(ISYMA,ISYM1)
        NI=NISH(ISYMI)
        ISYMIJ=Mul(ISYMI,ISYM12)
        NLST2=NLIST(ISYMI,ISYMIJ,15)
        IF(NA*NI*NLST2/=0) THEN
          LLST2=LLIST(ISYMI,ISYMIJ,15)
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
     &                X2(IX),DPT2(IDIT),Y(IY))
        END IF
      END DO
      END IF
C  -----------------------------------------------
      Case(15)

C  D&GP Two-el
      NLST1=NLIST(ISYM1,ISYM2, 8)
      NU=NASH(ISYM2)
      IF(NLST1*NU/=0) THEN
      LLST1=LLIST(ISYM1,ISYM2, 8)
      VAL1(1)= One
      VAL1(2)= One
      DO ISYMA=1,NSYM
        ISYMI=Mul(ISYMA,ISYM1)
        NI=NISH(ISYMI)
        ISYMAB=Mul(ISYMA,ISYM12)
        NLST2=NLIST(ISYMA,ISYMAB,16)
        IF(NI*NLST2/=0) THEN
          LLST2=LLIST(ISYMA,ISYMAB,16)
          VAL2(1)= SQRI2
          VAL2(2)= One
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
     &                X2(IX),DPT2(IDTA),Y(IY))
        END IF
      END DO
      END IF
C  -----------------------------------------------
      Case(16)

C  D&GM Two-el
      NLST1=NLIST(ISYM1,ISYM2, 8)
      NU=NASH(ISYM2)
      IF(NLST1*NU/=0) THEN
      LLST1=LLIST(ISYM1,ISYM2, 8)
      VAL1(1)=-One
      VAL1(2)= One

      DO ISYMA=1,NSYM
        ISYMI=Mul(ISYMA,ISYM1)
        NI=NISH(ISYMI)
        ISYMAB=Mul(ISYMA,ISYM12)
        NLST2=NLIST(ISYMA,ISYMAB,17)
        IF(NI*NLST2/=0) THEN
          LLST2=LLIST(ISYMA,ISYMAB,17)
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
     &                X2(IX),DPT2(IDTA),Y(IY))
        END IF
      END DO
      END IF
C  -----------------------------------------------
      Case(17)

C  EP&HP Two-el
      NA=NSSH(ISYM12)
      NLST1=NLIST(ISYM12,ISYM2,16)
      IF(NA*NLST1/=0) THEN
        LLST1=LLIST(ISYM12,ISYM2,16)
        VAL1(1)= SQRI2
        VAL1(2)= One
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
      END IF
C  -----------------------------------------------
      Case(18)

C  EM&HM Two-el
      NA=NSSH(ISYM12)
      NLST1=NLIST(ISYM12,ISYM2,17)
      IF(NA*NLST1/=0) THEN
        LLST1=LLIST(ISYM12,ISYM2,17)
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
      END IF
C  -----------------------------------------------
      Case(19)

C  FP&GP Two-el
      NI=NISH(ISYM12)
      NLST1=NLIST(ISYM1,ISYM2, 9)
      IF(NI*NLST1/=0) THEN
        LLST1=LLIST(ISYM1,ISYM2, 9)
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
      END IF
C  -----------------------------------------------
      Case(20)

C  FM&GM Two-el
      NI=NISH(ISYM12)
      NLST1=NLIST(ISYM1,ISYM2,10)
      IF(NI*NLST1/=0) THEN
        LLST1=LLIST(ISYM1,ISYM2,10)
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
      END IF
C  -----------------------------------------------
      Case(21)

C  GP&HP Two-el
      LLST1=LLIST(ISYM12,ISYM2,14)
      NLST1=NLIST(ISYM12,ISYM2,14)
      IF(NLST1/=0) THEN
        VAL1(1)= -SQRI2
        VAL1(2)= -One
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
      END IF
C  -----------------------------------------------
      Case(22)

C  GM&HM Two-el
      LLST1=LLIST(ISYM12,ISYM2,15)
      NLST1=NLIST(ISYM12,ISYM2,15)
      IF(NLST1/=0) THEN
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
      END IF
C  -----------------------------------------------
      Case(23)

C  D&HP One-el
      IF(ISYM1==1) THEN
      IOXIA=0
      DO ISYMI=1,NSYM
       NI=NISH(ISYMI)
       ISYMA=ISYMI
       NA=NSSH(ISYMA)
       ISYMJ=Mul(ISYMI,ISYM2)
       NJ=NISH(ISYMJ)
       ISYMB=ISYMJ
       NB=NSSH(ISYMB)
       NLST1=NLIST(ISYMI,ISYM2,14)
       NLST2=NLIST(ISYMA,ISYM2,16)
       IF(NI*NA*NJ*NB*NLST1*NLST2/=0) Then
         LLST1=LLIST(ISYMI,ISYM2,14)
         VAL1(1)=SQRI2
         VAL1(2)=One
         LLST2=LLIST(ISYMA,ISYM2,16)
         VAL2(1)=SQRI2
         VAL2(2)=One
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
     &                X1(IXIA),DPT2(IDJB),Y(IY))
       END IF
       IOXIA=IOXIA+NI*NA
      END DO
      END IF

C  ---------------------------
      Case(24)

C  D&HM One-el
      IF(ISYM1==1) THEN
      IOXIA=0
      DO ISYMI=1,NSYM
       NI=NISH(ISYMI)
       ISYMA=ISYMI
       NA=NSSH(ISYMA)
       ISYMJ=Mul(ISYMI,ISYM2)
       NJ=NISH(ISYMJ)
       ISYMB=ISYMJ
       NB=NSSH(ISYMB)
       NLST1=NLIST(ISYMI,ISYM2,15)
       NLST2=NLIST(ISYMA,ISYM2,17)
       IF(NI*NA*NJ*NB*NLST1*NLST2/=0) THEN
         LLST1=LLIST(ISYMI,ISYM2,15)
         VAL1(1)=SQR3*Half
         VAL1(2)=-VAL1(1)
         LLST2=LLIST(ISYMA,ISYM2,17)
         VAL2(1)= One
         VAL2(2)=-One
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
     &              X1(IXIA),DPT2(IDJB),Y(IY))
       END IF
       IOXIA=IOXIA+NI*NA
      END DO
      END IF
C  ---------------------------
      Case default
        WRITE(6,*)' INTERNAL ERROR: OffDns reached invalid KOD=',KOD
        Call Abend()
      End Select

      END SUBROUTINE OFFDNS
