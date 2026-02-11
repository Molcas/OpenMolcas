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
C SVC-20120305: This is a new parallel version of the SGM subroutine.
C Instead of the full array Y, distributed array indices lg_Y is passed,
C which refer to distributed arrays.

C Currently, only case H will be passed as a distributed array, since
C this is the largest array (about a factor of NI larger than case G)
C and it avoids communications of this array: if it is updated, we use
C only the local chunk. If is is used to update another case, that case
C is partially updated and then communicated.  The special routines are
C called MLTR1_GH, MLTR1_EH and MLTSCA_DH.

C The distributed arrays are (currently) layed out as 2-dimensional
C distributed arrays, with NAS rows and NIS columns (see subroutine
C RHS_ALLO in file par_rhs.f). The chunks are along the column indices,
C so each chunk has all the row indices (full columns).


      SUBROUTINE SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,
     &               X1,X2,lg_Y,LIST)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use Constants, only: One,  Half, Two, Three, Six
      use Fockof, only: IOFFIA, FIT, FTI, FIA, FAI, FTA, FAT
      use EQSOLV, only: nList, LList, IfCoup
      use Sigma_data, only: Val1, Val2, INCF1, INCF2, INCX1, INCX2,
     &                      INCX3, INCY1, INCY2, INCY3, LEN1, LEN2,
     &                      nLst1, nLst2
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: nIsh, nAsh, nSSh, nSym, nISUP,
     &                         nASUP, nIGEJ, nIGTJ, nAGEB, nAGTB, nTGEU,
     &                         nTUV, nTGTU, nTGEU
#ifdef _DEBUGPRINT_
      use caspt2_module, only: Cases
#endif

      IMPLICIT None

      integer(kind=iwp), intent(in):: IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2
      real(kind=wp), intent(inout) ::  X1(*), X2(*)
      integer(kind=iwp), intent(in) ::  LIST(*)

      integer(kind=iwp) IOFCD(8,8),IOFCEP(8,8),IOFCEM(8,8),IOFCGP(8,8),
     &                  IOFCGM(8,8)
C Various constants:
      real(kind=wp), parameter:: SQR2=SQRT(Two), SQR3=SQRT(Three),
     &                           SQR6=SQRT(Six), SQRI2=One/SQR2,
     &                           SQRI6=One/SQR6, SQR32=SQR3*SQRI2
      integer(kind=iwp) lg_Y, ICD, ICEM, ICEP, ICGM, ICGP, IJSYM, ISYM,
     &                  ISYM12, ISYMA, ISYMAB, ISYMB, ISYMI, ISYMIJ,
     &                  ISYMJ, IX, IXIA, IXTA, IXTI, IY, JSYM, JXOFF,
     &                  KOD, LLST1, LLST2, NA, NAS1, NAS2, NB, NFA, NFI,
     &                  NFT, NI, NIS1, NIS2, NJ, NT, NU

C Compute a contribution from a single block (ISYM2,ICASE2)
C of expansion vector to a single block (ISYM1,ICASE1)
C of a sigma vector, if IMLTOP=0. If IMLTOP=1, the role of the
C two blocks are reversed. It is assumed that ICASE2 is
C greater than ICASE1, hence the reverse option.

CFUE  IF( (ICASE1.EQ.1 .AND. ICASE2.EQ.2) .AND.
CFUE &    (ISYM1.EQ.1  .AND. ISYM2.EQ.1 )       ) IFTEST=1
CPAM      IF(ICASE1.EQ.5.AND.ICASE2.GT.11) IFTEST=1
#ifdef _DEBUGPRINT_
      WRITE(6,'(A,10I5)')' ENTERING SGM.'
      WRITE(6,'(A,10I5)')'       IMLTOP:',IMLTOP
      WRITE(6,'(A,10I5)')' ISYM1,ICASE1:',ISYM1,ICASE1
      WRITE(6,'(A,10I5)')' ISYM2,ICASE2:',ISYM2,ICASE2
      WRITE(6,'(A,A,A)') ' CASES: ',CASES(ICASE1),CASES(ICASE2)
#endif

      IF(IMLTOP.LT.0 .OR. IMLTOP.GT.2) THEN
        WRITE(6,*) 'Error in SGM: IMLTOP = ', IMLTOP
        CALL AbEnd()
      END IF

C SVC: IFCOUP is set in SIGMA_CASPT2
      KOD=IFCOUP(ICASE2,ICASE1)
      IF(KOD==0) RETURN

      ISYM12=Mul(ISYM1,ISYM2)
      NAS1=NASUP(ISYM1,ICASE1)
      NIS1=NISUP(ISYM1,ICASE1)
      NAS2=NASUP(ISYM2,ICASE2)
      NIS2=NISUP(ISYM2,ICASE2)
      DO ISYM=1,NSYM
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

C SVC: this is an extra check, since coupling cases that are 0 should
C not have entered the sgm subroutine

      SELECT CASE (KOD)
C  -----------------------------------------------
      CASE(1)
C ICASE1= 1
C ICASE2= 2

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
          INCF1=NISH(ISYM12)
          INCF2=1
          IY=1
          INCY1=1
          INCY2=NTGEU(ISYM2)
          CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                X1(IXTI),
     &                FIT(ISYM12)%A,
     &                GA_Arrays(lg_Y)%A(IY))
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
          INCF1=NISH(ISYM12)
          INCF2=1
          IY=1
          INCY1=1
          INCY2=NTGEU(ISYM2)
          CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                X2(IX),
     &                FIT(ISYM12)%A,
     &                GA_Arrays(lg_Y)%A(IY))
        END IF
C  -----------------------------------------------
      CASE (2)
C ICASE1= 1
C ICASE2= 3

C A&BM One-el
        NLST1=NLIST(ISYM1,ISYM2,13)
        NLST2=NLIST(ISYM1,ISYM2,15)
        IF(NLST1*NLST2/=0) THEN
          LLST1=LLIST(ISYM1,ISYM2,13)
          VAL1(1)= Three
          VAL1(2)=-Three
          LLST2=LLIST(ISYM1,ISYM2,15)
* Original:
*         VAL2(1)=-One
*         VAL2(2)= One
* Fix for sign error noted by Takeshi, May 2015:
          VAL2(1)= One
          VAL2(2)=-One
          IXTI=1
          INCX1=1
          INCX2=NASH(ISYM1)
          INCF1=NISH(ISYM12)
          INCF2=1
          IY=1
          INCY1=1
          INCY2=NTGTU(ISYM2)
          CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                X1(IXTI),
     &                FIT(ISYM12)%A,
     &                GA_Arrays(lg_Y)%A(IY))
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
*         VAL2(1)=-One
*         VAL2(2)= One
* Fix for sign error noted by Takeshi, May 2015:
          VAL2(1)= One
          VAL2(2)=-One
          IX=1
          INCX1=1
          INCX2=NTUV(ISYM1)
          INCF1=NISH(ISYM12)
          INCF2=1
          IY=1
          INCY1=1
          INCY2=NTGTU(ISYM2)
          CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                X2(IX),
     &                FIT(ISYM12)%A,
     &                GA_Arrays(lg_Y)%A(IY))
        END IF
C  -----------------------------------------------
      CASE (3)
C ICASE1= 1
C ICASE2= 5

C  A&D  Two-el
        NLST1=NLIST(ISYM1,ISYM2, 1)
        IF(NLST1/=0) THEN
          LLST1=LLIST(ISYM1,ISYM2, 1)
          VAL1(1)= One
          VAL1(2)= One
          IX=1
          INCX1=1
          INCX2=NAS1
          INCF1=NSSH(ISYM12)
          INCF2=1
          IY=1+NAS2*IOFCD(ISYM2,ISYM12)
          INCY1=1
          INCY2=NAS2
          INCY3=NAS2*NISH(ISYM1)
          LEN1=NISH(ISYM1)
          LEN2=NSSH(ISYM12)
          CALL MLTMV(IMLTOP,LIST(LLST1),
     &               X2(IX),
     &               FAT(ISYM12)%A,
     &               GA_Arrays(lg_Y)%A(IY))
        END IF
C  -----------------------------------------------
      CASE (4)
C ICASE1= 1
C ICASE2= 6

C  A&EP One-el
        IF(ISYM2.EQ.ISYM1) THEN
          NT=NASH(ISYM1)
          IF(NT/=0) THEN
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
              INCF1=NSSH(ISYMA)
              INCF2=1
              IY=1+NAS2*IOFCEP(ISYM2,ISYMA)
              INCY1=NT*NA
              INCY2=1
              INCY3=NT
              LEN1=NT
              LEN2=NA
              CALL MLTMV(IMLTOP,LIST(LLST1),
     &                   X1(IXTI),
     &                   FAI(ISYMA)%A,
     &                   GA_Arrays(lg_Y)%A(IY))
            END IF
            END DO
          END IF
        END IF
C  -----------------------------------------------
      CASE (5)
C ICASE1= 1
C ICASE2= 7

C  A&EM One-el
        IF(ISYM2.EQ.ISYM1) THEN
          NT=NASH(ISYM1)
          IF(NT/=0) THEN
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
              INCF1=NSSH(ISYMA)
              INCF2=1
              IY=1+NAS2*IOFCEM(ISYM2,ISYMA)
              INCY1=NT*NA
              INCY2=1
              INCY3=NT
              LEN1=NT
              LEN2=NA
              CALL MLTMV(IMLTOP,LIST(LLST1),
     &                   X1(IXTI),
     &                   FAI(ISYMA)%A,
     &                   GA_Arrays(lg_Y)%A(IY))
            END IF
            END DO
          END IF
        END IF
C  -----------------------------------------------
      CASE (6)
C ICASE1= 2
C ICASE2= 6

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
          INCF1=NSSH(ISYM12)
          INCF2=1
          IY=1+NAS2*IOFCEP(ISYM2,ISYM12)
          INCY1=1
          INCY2=NAS2*NA
          INCY3=NAS2
          LEN1=NIS1
          LEN2=NA
          CALL MLTMV(IMLTOP,LIST(LLST1),
     &               X2(IX),
     &               FAT(ISYM12)%A,
     &               GA_Arrays(lg_Y)%A(IY))
        END IF
C  -----------------------------------------------
      CASE (7)
C ICASE1= 3
C ICASE2= 7

C  BM&EM Two-el
        NA=NSSH(ISYM12)
        NLST1=NLIST(ISYM1,ISYM2,10)
        IF(NA*NLST1/=0) THEN
          LLST1=LLIST(ISYM1,ISYM2,10)
* Original:
*         VAL1(1)=-SQRI6
*         VAL1(2)= SQRI6
* Fix for sign error noted by Takeshi, May 2015:
          VAL1(1)= SQRI6
          VAL1(2)=-SQRI6
          IX=1
          INCX1=1
          INCX2=NAS1
          INCF1=NSSH(ISYM12)
          INCF2=1
          IY=1+NAS2*IOFCEM(ISYM2,ISYM12)
          INCY1=1
          INCY2=NAS2*NA
          INCY3=NAS2
          LEN1=NIS1
          LEN2=NA
          CALL MLTMV(IMLTOP,LIST(LLST1),
     &               X2(IX),
     &               FAT(ISYM12)%A,
     &               GA_Arrays(lg_Y)%A(IY))
        END IF
C  -----------------------------------------------
      CASE (8)
C ICASE1= 4
C ICASE2= 5

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
          INCF1=NI
          INCF2=1
          IY=1+NAS2*IOFCD(ISYM2,ISYM1)
          INCY1=1
          INCY2=NAS2*NI
          INCY3=NAS2
          LEN1=NSSH(ISYM1)
          LEN2=NI
          CALL MLTMV(IMLTOP,LIST(LLST1),
     &               X1(IXTA),
     &               FIT(ISYM12)%A,
     &               GA_Arrays(lg_Y)%A(IY))
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
          INCF1=NI
          INCF2=1
          IY=1+NAS2*IOFCD(ISYM2,ISYM1)
          INCY1=1
          INCY2=NAS2*NI
          INCY3=NAS2
          LEN1=NSSH(ISYM1)
          LEN2=NI
          CALL MLTMV(IMLTOP,LIST(LLST1),
     &               X2(IX),
     &               FIT(ISYM12)%A,
     &               GA_Arrays(lg_Y)%A(IY))
        END IF
C  -----------------------------------------------
      CASE (9)
C ICASE1= 4
C ICASE2= 8

C  C&FP One-el
        NLST1=NLIST(ISYM1,ISYM2,12)
        NLST2=NLIST(ISYM1,ISYM2,16)
        IF(NLST1*NLST2/=0) THEN
          LLST1=LLIST(ISYM1,ISYM2,12)
          VAL1(1)=-One
          VAL1(2)=-Two
          LLST2=LLIST(ISYM1,ISYM2,16)
          VAL2(1)=One
          VAL2(2)=SQR2
          IXTA=1
          INCX1=1
          INCX2=NASH(ISYM1)
          INCF1=1
          INCF2=NASH(ISYM12)
          IY=1
          INCY1=1
          INCY2=NTGEU(ISYM2)
          CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                X1(IXTA),
     &                FTA(ISYM12)%A,
     &                GA_Arrays(lg_Y)%A(IY))
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
          INCF1=1
          INCF2=NASH(ISYM12)
          IY=1
          INCY1=1
          INCY2=NTGEU(ISYM2)
          CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                X2(IX),
     &                FTA(ISYM12)%A,
     &                GA_Arrays(lg_Y)%A(IY))
        END IF
C  -----------------------------------------------
      CASE (10)
C ICASE1= 4
C ICASE2= 9

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
          INCF1=1
          INCF2=NASH(ISYM12)
          IY=1
          INCY1=1
          INCY2=NTGTU(ISYM2)
          CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                X1(IXTA),
     &                FTA(ISYM12)%A,
     &                GA_Arrays(lg_Y)%A(IY))
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
          INCF1=1
          INCF2=NASH(ISYM12)
          IY=1
          INCY1=1
          INCY2=NTGTU(ISYM2)
          CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                X2(IX),
     &                FTA(ISYM12)%A,
     &                GA_Arrays(lg_Y)%A(IY))
        END IF
C  -----------------------------------------------
      CASE (11)
C ICASE1= 4
C ICASE2= 10

C  C&GP One-el
        NT=NASH(ISYM1)
        IF(ISYM2==ISYM1 .AND.  NT/=0) THEN
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
            INCF1=NI
            INCF2=1
            IY=1+NAS2*IOFCGP(ISYM2,ISYMI)
            INCY1=NT*NI
            INCY2=1
            INCY3=NT
            LEN1=NT
            LEN2=NI
            CALL MLTMV(IMLTOP,LIST(LLST1),
     &                 X1(IXTA),
     &                 FIA(ISYMI)%A,
     &                 GA_Arrays(lg_Y)%A(IY))
          END IF
          END DO
        END IF
C  -----------------------------------------------
      CASE (12)
C ICASE1= 4
C ICASE2= 11

C  C&GM One-el
        NT=NASH(ISYM1)
        IF(ISYM2==ISYM1 .AND. NT/=0) THEN
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
            INCF1=NI
            INCF2=1
            IY=1+NAS2*IOFCGM(ISYM2,ISYMI)
            INCY1=NT*NI
            INCY2=1
            INCY3=NT
            LEN1=NT
            LEN2=NI
            CALL MLTMV(IMLTOP,LIST(LLST1),
     &                 X1(IXTA),
     &                 FIA(ISYMI)%A,
     &                 GA_Arrays(lg_Y)%A(IY))
          END IF
          END DO
        END IF
C  -----------------------------------------------
      CASE (13)
C ICASE1= 5
C ICASE2= 6

C  D&EP One-el
        NT=NASH(ISYM2)
        IF(ISYM1==1 .AND. NT/=0) THEN
          DO ISYMI=1,NSYM
          NI=NISH(ISYMI)
          NA=NSSH(ISYMI)
          ISYMIJ=Mul(ISYMI,ISYM2)
          NLST1=NLIST(ISYMI,ISYMIJ,14)
          IF(NI*NA*NLST1/=0) THEN
            LLST1=LLIST(ISYMI,ISYMIJ,14)
            VAL1(1)= SQRI2
            VAL1(2)= One
            IXIA=1+IOFFIA(ISYMI)
            INCX1=1
            INCX2=NI
            INCF1=NASH(ISYM2)
            INCF2=1
            IY=1+NAS2*IOFCEP(ISYM2,ISYMI)
            INCY1=NT*NA
            INCY2=NT
            INCY3=1
            LEN1=NA
            LEN2=NT
            CALL MLTMV(IMLTOP,LIST(LLST1),
     &                 X1(IXIA),
     &                 FTI(ISYM2)%A,
     &                 GA_Arrays(lg_Y)%A(IY))
          END IF
          END DO
        END IF

C  D&EP Two-el
        LLST1=LLIST(ISYM1,ISYM2, 7)
        NLST1=NLIST(ISYM1,ISYM2, 7)
        NU=NASH(ISYM2)
        IF(NLST1*NU/=0) THEN
          VAL1(1)=-One
          VAL1(2)=-One
          DO ISYMA=1,NSYM
          NA=NSSH(ISYMA)
          ISYMI=Mul(ISYMA,ISYM1)
          NI=NISH(ISYMI)
          ISYMIJ=Mul(ISYMI,ISYM12)
          NLST2=NLIST(ISYMI,ISYMIJ,14)
          IF(NA*NI*NLST2/=0) THEN
            LLST2=LLIST(ISYMI,ISYMIJ,14)
            VAL2(1)= SQRI2
            VAL2(2)= One
            IX=1+NAS1*IOFCD(ISYM1,ISYMA)
            INCX1=1
            INCX2=NAS1
            INCX3=NAS1*NI
            INCF1=NISH(ISYM12)
            INCF2=1
            IY=1+NU*IOFCEP(ISYM2,ISYMA)
            INCY1=1
            INCY2=NU*NA
            INCY3=NU
            LEN1=NA
            CALL MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                  X2(IX),
     &                  FIT(ISYM12)%A,
     &                  GA_Arrays(lg_Y)%A(IY))
          END IF
          END DO
        END IF
C  -----------------------------------------------
      CASE (14)
C ICASE1= 5
C ICASE2= 7

C  D&EM One-el
        NT=NASH(ISYM2)
        IF(ISYM1==1 .AND. NT/=0) THEN
          DO ISYMI=1,NSYM
          NI=NISH(ISYMI)
          NA=NSSH(ISYMI)
          ISYMIJ=Mul(ISYMI,ISYM2)
          NLST1=NLIST(ISYMI,ISYMIJ,15)
          IF(NI*NA*NLST1/=0) THEN
            LLST1=LLIST(ISYMI,ISYMIJ,15)
            VAL1(1)= SQR32
            VAL1(2)=-SQR32
            IXIA=1+IOFFIA(ISYMI)
            INCX1=1
            INCX2=NI
            INCF1=NASH(ISYM2)
            INCF2=1
            IY=1+NAS2*IOFCEM(ISYM2,ISYMI)
            INCY1=NT*NA
            INCY2=NT
            INCY3=1
            LEN1=NA
            LEN2=NT
            CALL MLTMV(IMLTOP,LIST(LLST1),
     &                 X1(IXIA),
     &                 FTI(ISYM2)%A,
     &                 GA_Arrays(lg_Y)%A(IY))
          END IF
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
              INCF1=NISH(ISYM12)
              INCF2=1
              IY=1+NU*IOFCEM(ISYM2,ISYMA)
              INCY1=1
              INCY2=NU*NA
              INCY3=NU
              LEN1=NA
              CALL MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                          X2(IX),
     &                          FIT(ISYM12)%A,
     &                          GA_Arrays(lg_Y)%A(IY))
            END IF
            END DO
        END IF
C  -----------------------------------------------
      CASE (15)
C ICASE1= 5
C ICASE2= 10

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
              INCF1=1
              INCF2=NASH(ISYM12)
              IY=1+NU*IOFCGP(ISYM2,ISYMI)
              INCY1=1
              INCY2=NU*NI
              INCY3=NU
              LEN1=NI
              CALL MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                    X2(IX),
     &                    FTA(ISYM12)%A,
     &                    GA_Arrays(lg_Y)%A(IY))
            END IF
            END DO
        END IF
C  -----------------------------------------------
      CASE (16)
C ICASE1= 5
C ICASE2= 11

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
              INCF1=1
              INCF2=NASH(ISYM12)
              IY=1+NU*IOFCGM(ISYM2,ISYMI)
              INCY1=1
              INCY2=NU*NI
              INCY3=NU
              LEN1=NI
              CALL MLTDXP(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                    X2(IX),
     &                    FTA(ISYM12)%A,
     &                    GA_Arrays(lg_Y)%A(IY))
            END IF
            END DO
        END IF
C  -----------------------------------------------
      CASE (17)
C ICASE1= 6
C ICASE2= 12

C  EP&HP Two-el
        NA=NSSH(ISYM12)
        NLST1=NLIST(ISYM12,ISYM2,16)
        IF(NA*NLST1/=0) THEN
          LLST1=LLIST(ISYM12,ISYM2,16)
          VAL1(1)= SQRI2
          VAL1(2)= One
          JXOFF=IOFCEP(ISYM1,ISYM12)
          INCX3=NAS1*NA
          NFT=NASH(ISYM1)
          NFA=NSSH(ISYM1)
          CALL PMLTR1(KOD,IMLTOP,LIST(LLST1),
     &                X2,NAS1,NIS1,JXOFF,
     &                FTA(ISYM1)%A,NFT,NFA,
     &                lg_Y,NAS2,NIS2)
        END IF
C  -----------------------------------------------
      CASE (18)
C ICASE1= 7
C ICASE2= 13

C  EM&HM Two-el
        NA=NSSH(ISYM12)
        NLST1=NLIST(ISYM12,ISYM2,17)
        IF(NA*NLST1/=0) THEN
          LLST1=LLIST(ISYM12,ISYM2,17)
          VAL1(1)= SQRI2
          VAL1(2)=-SQRI2
          JXOFF=IOFCEM(ISYM1,ISYM12)
          INCX3=NAS1*NA
          NFT=NASH(ISYM1)
          NFA=NSSH(ISYM1)
          CALL PMLTR1(KOD,IMLTOP,LIST(LLST1),
     &                X2,NAS1,NIS1,JXOFF,
     &                FTA(ISYM1)%A,NFT,NFA,
     &                lg_Y,NAS2,NIS2)
        END IF
C  -----------------------------------------------
      CASE (19)
C ICASE1= 8
C ICASE2= 10

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
          INCF1=NI
          INCF2=1
          IY=1+NAS2*IOFCGP(ISYM2,ISYM12)
          INCY1=1
          INCY2=NAS2*NI
          INCY3=NAS2
          LEN1=NIS1
          LEN2=NI
          CALL MLTMV(IMLTOP,LIST(LLST1),
     &               X2(IX),
     &               FIT(ISYM12)%A,
     &               GA_Arrays(lg_Y)%A(IY))
        END IF
C  -----------------------------------------------
      CASE (20)
C ICASE1= 9
C ICASE2= 11

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
          INCF1=NI
          INCF2=1
          IY=1+NAS2*IOFCGM(ISYM2,ISYM12)
          INCY1=1
          INCY2=NAS2*NI
          INCY3=NAS2
          LEN1=NIS1
          LEN2=NI
          CALL MLTMV(IMLTOP,LIST(LLST1),
     &               X2(IX),
     &               FIT(ISYM12)%A,
     &               GA_Arrays(lg_Y)%A(IY))
        END IF
C  -----------------------------------------------
      CASE (21)
C ICASE1= 10
C ICASE2= 12

C  GP&HP Two-el
        NLST1=NLIST(ISYM12,ISYM2,14)
        IF(NLST1/=0) THEN
          LLST1=LLIST(ISYM12,ISYM2,14)
          VAL1(1)= -SQRI2
          VAL1(2)= -One
          JXOFF=IOFCGP(ISYM1,ISYM12)
          INCX3=NAS1*NISH(ISYM12)
          NFT=NASH(ISYM1)
          NFI=NISH(ISYM1)
          CALL PMLTR1(KOD,IMLTOP,LIST(LLST1),
     &                X2,NAS1,NIS1,JXOFF,
     &                FTI(ISYM1)%A,NFT,NFI,
     &                lg_Y,NAS2,NIS2)
        END IF
C  -----------------------------------------------
      CASE (22)
C ICASE1= 11
C ICASE2= 13

C  GM&HM Two-el
        NLST1=NLIST(ISYM12,ISYM2,15)
        IF(NLST1/=0) THEN
          LLST1=LLIST(ISYM12,ISYM2,15)
          VAL1(1)= SQRI2
          VAL1(2)=-SQRI2
          JXOFF=IOFCGM(ISYM1,ISYM12)
          INCX3=NAS1*NISH(ISYM12)
          NFT=NASH(ISYM1)
          NFI=NISH(ISYM1)
          CALL PMLTR1(KOD,IMLTOP,LIST(LLST1),
     &                X2,NAS1,NIS1,JXOFF,
     &                FTI(ISYM1)%A,NFT,NFI,
     &                lg_Y,NAS2,NIS2)
        END IF
C  -----------------------------------------------
      CASE (23)
C ICASE1= 5
C ICASE2= 12

C  D&HP One-el
        IF(ISYM1==1) THEN
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
          IF(NI*NA*NJ*NB*NLST1*NLST2/=0) THEN
            LLST1=LLIST(ISYMI,ISYM2,14)
            VAL1(1)=SQRI2
            VAL1(2)=One
            LLST2=LLIST(ISYMA,ISYM2,16)
            VAL2(1)=SQRI2
            VAL2(2)=One
            IXIA=1+IOFFIA(ISYMI)
            INCX1=1
            INCX2=NI
            CALL PMLTSCA(KOD,IMLTOP,
     &                   LIST(LLST1),LIST(LLST2),
     &                   X1(IXIA),NI,NA,
     &                   FIA(ISYMB)%A,NJ,NB,
     &                   lg_Y,NAS2,NIS2)
          END IF
          END DO
        END IF
C ---------------------------
      CASE (24)
C ICASE1= 5
C ICASE2= 13

C  D&HM One-el
        IF(ISYM1==1) THEN
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
            IXIA=1+IOFFIA(ISYMI)
            INCX1=1
            INCX2=NI
            CALL PMLTSCA(KOD,IMLTOP,
     &                   LIST(LLST1),LIST(LLST2),
     &                   X1(IXIA),NI,NA,
     &                   FIA(ISYMB)%A,NJ,NB,
     &                   lg_Y,NAS2,NIS2)
          END IF
          END DO
        END IF
C ---------------------------
      CASE DEFAULT
        WRITE(6,*)' INTERNAL ERROR: SGM reached invalid KOD=',KOD
        Call Abend()
      END SELECT

      END SUBROUTINE SGM
