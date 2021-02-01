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
C Instead of the full arrays X2 and Y, distributed array indices lg_X
C and lg_Y are passed, they can be either indices into the WORK array or
C they can refer to distributed arrays.

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
     &               X1,lg_X,lg_Y,LIST)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
      DIMENSION X1(*)
      DIMENSION LIST(*)
      DIMENSION IOFCD(8,8),IOFCEP(8,8),IOFCEM(8,8),IOFCGP(8,8),
     &          IOFCGM(8,8)
#include "sigma.fh"
      COMMON /CPLCAS/ IFCOUP(MXCASE,MXCASE)
      COMMON /FOCKOF/ LFIT,LFIA,LFTA,LFTI,LFAI,LFAT,
     &                IOFFIT(8),IOFFIA(8),IOFFTA(8)
C Various constants:
      SQR2=SQRT(2.0D00)
      SQR3=SQRT(3.0D00)
      SQR6=SQRT(6.0D00)
      SQRI2=1.0D00/SQR2
      SQRI6=1.0D00/SQR6
      SQR32=SQR3*SQRI2

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

      ISYM12=MUL(ISYM1,ISYM2)
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

C SVC: this is an extra check, since coupling cases that are 0 should
C not have entered the sgm subroutine
      IF(KOD.EQ.0) THEN
        RETURN
C  -----------------------------------------------
      ELSE IF (KOD.EQ.1) THEN
C ICASE1= 1
C ICASE2= 2

C  A&BP One-el
CFUE  Call GetMem('CX','Check',' ',iDummy,iDummy)
CTEST      WRITE(*,*) ' A&BP One-el'
        LLST1=LLIST(ISYM1,ISYM2,12)
        NLST1=NLIST(ISYM1,ISYM2,12)
        IF(NLST1.NE.0) THEN
          VAL1(1)=1.0D00
          VAL1(2)=2.0D00
          LLST2=LLIST(ISYM1,ISYM2,14)
          NLST2=NLIST(ISYM1,ISYM2,14)
          IF(NLST2.NE.0) THEN
            VAL2(1)=1.0D00
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
     &                  X1(IXTI),
     &                  WORK(LFIT+IOFFIT(ISYM12)),
     &                  WORK(lg_Y+IY-1))
          END IF
        END IF

C  A&BP Two-el
        LLST1=LLIST(ISYM1,ISYM2, 3)
        NLST1=NLIST(ISYM1,ISYM2, 3)
        IF(NLST1.NE.0) THEN
          VAL1(1)=-1.0D00
          VAL1(2)=-2.0D00
          LLST2=LLIST(ISYM1,ISYM2,14)
          NLST2=NLIST(ISYM1,ISYM2,14)
          IF(NLST2.NE.0) THEN
            VAL2(1)=1.0D00
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
     &                  WORK(lg_X+IX-1),
     &                  WORK(LFIT+IOFFIT(ISYM12)),
     &                  WORK(lg_Y+IY-1))
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.2) THEN
C ICASE1= 1
C ICASE2= 3

C A&BM One-el
        LLST1=LLIST(ISYM1,ISYM2,13)
        NLST1=NLIST(ISYM1,ISYM2,13)
        IF(NLST1.NE.0) THEN
          VAL1(1)= 3.0D00
          VAL1(2)=-3.0D00
          LLST2=LLIST(ISYM1,ISYM2,15)
          NLST2=NLIST(ISYM1,ISYM2,15)
          IF(NLST2.NE.0) THEN
* Original:
*           VAL2(1)=-1.0D00
*           VAL2(2)= 1.0D00
* Fix for sign error noted by Takeshi, May 2015:
            VAL2(1)= 1.0D00
            VAL2(2)=-1.0D00
            IXTI=1
            INCX1=1
            INCX2=NASH(ISYM1)
            INCF1=NISH(ISYM12)
            INCF2=1
            IY=1
            INCY1=1
            INCY2=NTGTU(ISYM2)
            CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                  X1(IXTI),
     &                  WORK(LFIT+IOFFIT(ISYM12)),
     &                  WORK(lg_Y+IY-1))
          END IF
        END IF

C A&BM Two-el
        LLST1=LLIST(ISYM1,ISYM2, 4)
        NLST1=NLIST(ISYM1,ISYM2, 4)
        IF(NLST1.NE.0) THEN
          VAL1(1)=-1.0D00
          VAL1(2)= 1.0D00
          LLST2=LLIST(ISYM1,ISYM2,15)
          NLST2=NLIST(ISYM1,ISYM2,15)
          IF(NLST2.NE.0) THEN
* Original:
*           VAL2(1)=-1.0D00
*           VAL2(2)= 1.0D00
* Fix for sign error noted by Takeshi, May 2015:
            VAL2(1)= 1.0D00
            VAL2(2)=-1.0D00
            IX=1
            INCX1=1
            INCX2=NTUV(ISYM1)
            INCF1=NISH(ISYM12)
            INCF2=1
            IY=1
            INCY1=1
            INCY2=NTGTU(ISYM2)
            CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                  WORK(lg_X+IX-1),
     &                  WORK(LFIT+IOFFIT(ISYM12)),
     &                  WORK(lg_Y+IY-1))
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.3) THEN
C ICASE1= 1
C ICASE2= 5

C  A&D  Two-el
        LLST1=LLIST(ISYM1,ISYM2, 1)
        NLST1=NLIST(ISYM1,ISYM2, 1)
        IF(NLST1.NE.0) THEN
          VAL1(1)= 1.0D00
          VAL1(2)= 1.0D00
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
     &               WORK(lg_X+IX-1),
     &               WORK(LFAT+IOFFTA(ISYM12)),
     &               WORK(lg_Y+IY-1))
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.4) THEN
C ICASE1= 1
C ICASE2= 6

C  A&EP One-el
        IF(ISYM2.EQ.ISYM1) THEN
          NT=NASH(ISYM1)
          IF(NT.NE.0) THEN
            DO ISYMIJ=1,NSYM
            ISYMA=MUL(ISYMIJ,ISYM1)
            NA=NSSH(ISYMA)
            IF(NA.NE.0) THEN
              LLST1=LLIST(ISYM1,ISYMIJ,14)
              NLST1=NLIST(ISYM1,ISYMIJ,14)
              IF(NLST1.NE.0) THEN
                VAL1(1)= 1.0D00
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
     &                     X1(IXTI),
     &                     WORK(LFAI+IOFFIA(ISYMA)),
     &                     WORK(lg_Y+IY-1))
              END IF
            END IF
            END DO
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.5) THEN
C ICASE1= 1
C ICASE2= 7

C  A&EM One-el
        IF(ISYM2.EQ.ISYM1) THEN
          NT=NASH(ISYM1)
          IF(NT.NE.0) THEN
            DO ISYMIJ=1,NSYM
            ISYMA=MUL(ISYMIJ,ISYM1)
            NA=NSSH(ISYMA)
            IF(NA.NE.0) THEN
              LLST1=LLIST(ISYM1,ISYMIJ,15)
              NLST1=NLIST(ISYM1,ISYMIJ,15)
              IF(NLST1.NE.0) THEN
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
     &                     X1(IXTI),
     &                     WORK(LFAI+IOFFIA(ISYMA)),
     &                     WORK(lg_Y+IY-1))
              END IF
            END IF
            END DO
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.6) THEN
C ICASE1= 2
C ICASE2= 6

C  BP&EP Two-el
        NA=NSSH(ISYM12)
        IF(NA.NE.0) THEN
          LLST1=LLIST(ISYM1,ISYM2, 9)
          NLST1=NLIST(ISYM1,ISYM2, 9)
          IF(NLST1.NE.0) THEN
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
     &                 WORK(lg_X+IX-1),
     &                 WORK(LFAT+IOFFTA(ISYM12)),
     &                 WORK(lg_Y+IY-1))
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.7) THEN
C ICASE1= 3
C ICASE2= 7

C  BM&EM Two-el
        NA=NSSH(ISYM12)
        IF(NA.NE.0) THEN
          LLST1=LLIST(ISYM1,ISYM2,10)
          NLST1=NLIST(ISYM1,ISYM2,10)
          IF(NLST1.NE.0) THEN
* Original:
*           VAL1(1)=-SQRI6
*           VAL1(2)= SQRI6
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
     &                 WORK(lg_X+IX-1),
     &                 WORK(LFAT+IOFFTA(ISYM12)),
     &                 WORK(lg_Y+IY-1))
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.8) THEN
C ICASE1= 4
C ICASE2= 5

C  C&D  One-el
        NI=NISH(ISYM12)
        IF(NI.NE.0) THEN
          LLST1=LLIST(ISYM1,ISYM2,11)
          NLST1=NLIST(ISYM1,ISYM2,11)
          IF(NLST1.NE.0) THEN
            VAL1(1)= 2.0D00
            VAL1(2)= 1.0D00
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
     &                 X1(IXTA),
     &                 WORK(LFIT+IOFFIT(ISYM12)),
     &                 WORK(lg_Y+IY-1))
          END IF
        END IF

C  C&D  Two-el
        NI=NISH(ISYM12)
        IF(NI.NE.0) THEN
          LLST1=LLIST(ISYM1,ISYM2, 2)
          NLST1=NLIST(ISYM1,ISYM2, 2)
          IF(NLST1.NE.0) THEN
            VAL1(1)=-1.0D00
            VAL1(2)=-1.0D00
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
     &                 WORK(lg_X+IX-1),
     &                 WORK(LFIT+IOFFIT(ISYM12)),
     &                 WORK(lg_Y+IY-1))
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.9) THEN
C ICASE1= 4
C ICASE2= 8

C  C&FP One-el
        LLST1=LLIST(ISYM1,ISYM2,12)
        NLST1=NLIST(ISYM1,ISYM2,12)
        IF(NLST1.NE.0) THEN
          VAL1(1)=-1.0D00
          VAL1(2)=-2.0D00
          LLST2=LLIST(ISYM1,ISYM2,16)
          NLST2=NLIST(ISYM1,ISYM2,16)
          IF(NLST2.NE.0) THEN
            VAL2(1)=1.0D00
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
     &                  X1(IXTA),
     &                  WORK(LFTA+IOFFTA(ISYM12)),
     &                  WORK(lg_Y+IY-1))
          END IF
        END IF

C  C&FP Two-el
        LLST1=LLIST(ISYM1,ISYM2, 5)
        NLST1=NLIST(ISYM1,ISYM2, 5)
        IF(NLST1.NE.0) THEN
          VAL1(1)= 1.0D00
          VAL1(2)= 2.0D00
          LLST2=LLIST(ISYM1,ISYM2,16)
          NLST2=NLIST(ISYM1,ISYM2,16)
          IF(NLST2.NE.0) THEN
            VAL2(1)=1.0D00
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
     &                  WORK(lg_X+IX-1),
     &                  WORK(LFTA+IOFFTA(ISYM12)),
     &                  WORK(lg_Y+IY-1))
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.10) THEN
C ICASE1= 4
C ICASE2= 9

C  C&FM One-el
        LLST1=LLIST(ISYM1,ISYM2,13)
        NLST1=NLIST(ISYM1,ISYM2,13)
        IF(NLST1.NE.0) THEN
          VAL1(1)=-1.0D00
          VAL1(2)= 1.0D00
          LLST2=LLIST(ISYM1,ISYM2,17)
          NLST2=NLIST(ISYM1,ISYM2,17)
          IF(NLST2.NE.0) THEN
            VAL2(1)=1.0D00
            VAL2(2)=-1.0D00
            IXTA=1
            INCX1=1
            INCX2=NASH(ISYM1)
            INCF1=1
            INCF2=NASH(ISYM12)
            IY=1
            INCY1=1
            INCY2=NTGTU(ISYM2)
            CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                  X1(IXTA),
     &                  WORK(LFTA+IOFFTA(ISYM12)),
     &                  WORK(lg_Y+IY-1))
          END IF
        END IF

C  C&FM Two-el
        LLST1=LLIST(ISYM1,ISYM2, 6)
        NLST1=NLIST(ISYM1,ISYM2, 6)
        IF(NLST1.NE.0) THEN
          VAL1(1)=-1.0D00
          VAL1(2)= 1.0D00
          LLST2=LLIST(ISYM1,ISYM2,17)
          NLST2=NLIST(ISYM1,ISYM2,17)
          IF(NLST2.NE.0) THEN
            VAL2(1)=1.0D00
            VAL2(2)=-1.0D00
            IX=1
            INCX1=1
            INCX2=NTUV(ISYM1)
            INCF1=1
            INCF2=NASH(ISYM12)
            IY=1
            INCY1=1
            INCY2=NTGTU(ISYM2)
            CALL MLTSCA(IMLTOP,LIST(LLST1),LIST(LLST2),
     &                  WORK(lg_X+IX-1),
     &                  WORK(LFTA+IOFFTA(ISYM12)),
     &                  WORK(lg_Y+IY-1))
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.11) THEN
C ICASE1= 4
C ICASE2= 10

C  C&GP One-el
        IF(ISYM2.EQ.ISYM1) THEN
          NT=NASH(ISYM1)
          IF(NT.NE.0) THEN
            DO ISYMAB=1,NSYM
            ISYMI=MUL(ISYMAB,ISYM1)
            NI=NISH(ISYMI)
            IF(NI.NE.0) THEN
              LLST1=LLIST(ISYM1,ISYMAB,16)
              NLST1=NLIST(ISYM1,ISYMAB,16)
              IF(NLST1.NE.0) THEN
                VAL1(1)= SQRI2
                VAL1(2)= 1.0D00
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
     &                     X1(IXTA),
     &                     WORK(LFIA+IOFFIA(ISYMI)),
     &                     WORK(lg_Y+IY-1))
              END IF
            END IF
            END DO
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.12) THEN
C ICASE1= 4
C ICASE2= 11

C  C&GM One-el
        IF(ISYM2.EQ.ISYM1) THEN
          NT=NASH(ISYM1)
          IF(NT.NE.0) THEN
            DO ISYMAB=1,NSYM
            ISYMI=MUL(ISYMAB,ISYM1)
            NI=NISH(ISYMI)
            IF(NI.NE.0) THEN
              LLST1=LLIST(ISYM1,ISYMAB,17)
              NLST1=NLIST(ISYM1,ISYMAB,17)
              IF(NLST1.NE.0) THEN
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
     &                     X1(IXTA),
     &                     WORK(LFIA+IOFFIA(ISYMI)),
     &                     WORK(lg_Y+IY-1))
              END IF
            END IF
            END DO
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.13) THEN
C ICASE1= 5
C ICASE2= 6

C  D&EP One-el
        IF(ISYM1.EQ.1) THEN
          NT=NASH(ISYM2)
          IF(NT.NE.0) THEN
            DO ISYMI=1,NSYM
            NI=NISH(ISYMI)
            IF(NI.NE.0) THEN
              NA=NSSH(ISYMI)
              IF(NA.NE.0) THEN
                ISYMIJ=MUL(ISYMI,ISYM2)
                LLST1=LLIST(ISYMI,ISYMIJ,14)
                NLST1=NLIST(ISYMI,ISYMIJ,14)
                IF(NLST1.NE.0) THEN
                  VAL1(1)= SQRI2
                  VAL1(2)= 1.0D00
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
     &                       X1(IXIA),
     &                       WORK(LFTI+IOFFIT(ISYM2)),
     &                       WORK(lg_Y+IY-1))
                END IF
              END IF
            END IF
            END DO
          END IF
        END IF

C  D&EP Two-el
        LLST1=LLIST(ISYM1,ISYM2, 7)
        NLST1=NLIST(ISYM1,ISYM2, 7)
        IF(NLST1.NE.0) THEN
          VAL1(1)=-1.0D00
          VAL1(2)=-1.0D00
          NU=NASH(ISYM2)
          IF(NU.NE.0) THEN
            DO ISYMA=1,NSYM
            NA=NSSH(ISYMA)
            IF(NA.NE.0) THEN
              ISYMI=MUL(ISYMA,ISYM1)
              NI=NISH(ISYMI)
              IF(NI.NE.0) THEN
                ISYMIJ=MUL(ISYMI,ISYM12)
                LLST2=LLIST(ISYMI,ISYMIJ,14)
                NLST2=NLIST(ISYMI,ISYMIJ,14)
                IF(NLST2.NE.0) THEN
                  VAL2(1)= SQRI2
                  VAL2(2)= 1.0D00
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
     &                        WORK(lg_X+IX-1),
     &                        WORK(LFIT+IOFFIT(ISYM12)),
     &                        WORK(lg_Y+IY-1))
                END IF
              END IF
            END IF
            END DO
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.14) THEN
C ICASE1= 5
C ICASE2= 7

C  D&EM One-el
        IF(ISYM1.EQ.1) THEN
          NT=NASH(ISYM2)
          IF(NT.NE.0) THEN
            DO ISYMI=1,NSYM
            NI=NISH(ISYMI)
            IF(NI.NE.0) THEN
              NA=NSSH(ISYMI)
              IF(NI.NE.0) THEN
                ISYMIJ=MUL(ISYMI,ISYM2)
                LLST1=LLIST(ISYMI,ISYMIJ,15)
                NLST1=NLIST(ISYMI,ISYMIJ,15)
                IF(NLST1.NE.0) THEN
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
     &                       X1(IXIA),
     &                       WORK(LFTI+IOFFIT(ISYM2)),
     &                       WORK(lg_Y+IY-1))
                END IF
              END IF
            END IF
            END DO
          END IF
        END IF

C  D&EM Two-el
        LLST1=LLIST(ISYM1,ISYM2, 7)
        NLST1=NLIST(ISYM1,ISYM2, 7)
        IF(NLST1.NE.0) THEN
          VAL1(1)=-1.0D00
          VAL1(2)= 1.0D00
          NU=NASH(ISYM2)
          IF(NU.NE.0) THEN
            DO ISYMA=1,NSYM
            NA=NSSH(ISYMA)
            IF(NA.NE.0) THEN
              ISYMI=MUL(ISYMA,ISYM1)
              NI=NISH(ISYMI)
              IF(NI.NE.0) THEN
                ISYMIJ=MUL(ISYMI,ISYM12)
                LLST2=LLIST(ISYMI,ISYMIJ,15)
                NLST2=NLIST(ISYMI,ISYMIJ,15)
                IF(NLST2.NE.0) THEN
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
     &                              WORK(lg_X+IX-1),
     &                              WORK(LFIT+IOFFIT(ISYM12)),
     &                              WORK(lg_Y+IY-1))
                END IF
              END IF
            END IF
            END DO
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.15) THEN
C ICASE1= 5
C ICASE2= 10

C  D&GP Two-el
        LLST1=LLIST(ISYM1,ISYM2, 8)
        NLST1=NLIST(ISYM1,ISYM2, 8)
        IF(NLST1.NE.0) THEN
          VAL1(1)= 1.0D00
          VAL1(2)= 1.0D00
          NU=NASH(ISYM2)
          IF(NU.NE.0) THEN
            DO ISYMA=1,NSYM
            ISYMI=MUL(ISYMA,ISYM1)
            NI=NISH(ISYMI)
            IF(NI.NE.0) THEN
              ISYMAB=MUL(ISYMA,ISYM12)
              LLST2=LLIST(ISYMA,ISYMAB,16)
              NLST2=NLIST(ISYMA,ISYMAB,16)
              IF(NLST2.NE.0) THEN
                VAL2(1)= SQRI2
                VAL2(2)= 1.0D00
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
     &                      WORK(lg_X+IX-1),
     &                      WORK(LFTA+IOFFTA(ISYM12)),
     &                      WORK(lg_Y+IY-1))
              END IF
            END IF
            END DO
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.16) THEN
C ICASE1= 5
C ICASE2= 11

C  D&GM Two-el
        LLST1=LLIST(ISYM1,ISYM2, 8)
        NLST1=NLIST(ISYM1,ISYM2, 8)
        IF(NLST1.NE.0) THEN
          VAL1(1)=-1.0D00
          VAL1(2)= 1.0D00
          NU=NASH(ISYM2)
          IF(NU.NE.0) THEN
            DO ISYMA=1,NSYM
            ISYMI=MUL(ISYMA,ISYM1)
            NI=NISH(ISYMI)
            IF(NI.NE.0) THEN
              ISYMAB=MUL(ISYMA,ISYM12)
              LLST2=LLIST(ISYMA,ISYMAB,17)
              NLST2=NLIST(ISYMA,ISYMAB,17)
              IF(NLST2.NE.0) THEN
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
     &                      WORK(lg_X+IX-1),
     &                      WORK(LFTA+IOFFTA(ISYM12)),
     &                      WORK(lg_Y+IY-1))
              END IF
            END IF
            END DO
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.17) THEN
C ICASE1= 6
C ICASE2= 12

C  EP&HP Two-el
        NA=NSSH(ISYM12)
        IF(NA.NE.0) THEN
          LLST1=LLIST(ISYM12,ISYM2,16)
          NLST1=NLIST(ISYM12,ISYM2,16)
          IF(NLST1.NE.0) THEN
            VAL1(1)= SQRI2
            VAL1(2)= 1.0D00
            JXOFF=IOFCEP(ISYM1,ISYM12)
            INCX3=NAS1*NA
            NFT=NASH(ISYM1)
            NFA=NSSH(ISYM1)
            CALL PMLTR1(KOD,IMLTOP,LIST(LLST1),
     &                  lg_X,NAS1,NIS1,JXOFF,
     &                  WORK(LFTA+IOFFTA(ISYM1)),NFT,NFA,
     &                  lg_Y,NAS2,NIS2)
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.18) THEN
C ICASE1= 7
C ICASE2= 13

C  EM&HM Two-el
        NA=NSSH(ISYM12)
        IF(NA.NE.0) THEN
          LLST1=LLIST(ISYM12,ISYM2,17)
          NLST1=NLIST(ISYM12,ISYM2,17)
          IF(NLST1.NE.0) THEN
            VAL1(1)= SQRI2
            VAL1(2)=-SQRI2
            JXOFF=IOFCEM(ISYM1,ISYM12)
            INCX3=NAS1*NA
            NFT=NASH(ISYM1)
            NFA=NSSH(ISYM1)
            CALL PMLTR1(KOD,IMLTOP,LIST(LLST1),
     &                  lg_X,NAS1,NIS1,JXOFF,
     &                  WORK(LFTA+IOFFTA(ISYM1)),NFT,NFA,
     &                  lg_Y,NAS2,NIS2)
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.19) THEN
C ICASE1= 8
C ICASE2= 10

C  FP&GP Two-el
        NI=NISH(ISYM12)
        IF(NI.NE.0) THEN
          LLST1=LLIST(ISYM1,ISYM2, 9)
          NLST1=NLIST(ISYM1,ISYM2, 9)
          IF(NLST1.NE.0) THEN
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
     &                 WORK(lg_X+IX-1),
     &                 WORK(LFIT+IOFFIT(ISYM12)),
     &                 WORK(lg_Y+IY-1))
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.20) THEN
C ICASE1= 9
C ICASE2= 11

C  FM&GM Two-el
        NI=NISH(ISYM12)
        IF(NI.NE.0) THEN
          LLST1=LLIST(ISYM1,ISYM2,10)
          NLST1=NLIST(ISYM1,ISYM2,10)
          IF(NLST1.NE.0) THEN
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
     &                 WORK(lg_X+IX-1),
     &                 WORK(LFIT+IOFFIT(ISYM12)),
     &                 WORK(lg_Y+IY-1))
          END IF
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.21) THEN
C ICASE1= 10
C ICASE2= 12

C  GP&HP Two-el
        LLST1=LLIST(ISYM12,ISYM2,14)
        NLST1=NLIST(ISYM12,ISYM2,14)
        IF(NLST1.NE.0) THEN
          VAL1(1)= -SQRI2
          VAL1(2)= -1.0D00
          JXOFF=IOFCGP(ISYM1,ISYM12)
          INCX3=NAS1*NISH(ISYM12)
          NFT=NASH(ISYM1)
          NFI=NISH(ISYM1)
          CALL PMLTR1(KOD,IMLTOP,LIST(LLST1),
     &                lg_X,NAS1,NIS1,JXOFF,
     &                WORK(LFTI+IOFFIT(ISYM1)),NFT,NFI,
     &                lg_Y,NAS2,NIS2)
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.22) THEN
C ICASE1= 11
C ICASE2= 13

C  GM&HM Two-el
        LLST1=LLIST(ISYM12,ISYM2,15)
        NLST1=NLIST(ISYM12,ISYM2,15)
        IF(NLST1.NE.0) THEN
          VAL1(1)= SQRI2
          VAL1(2)=-SQRI2
          JXOFF=IOFCGM(ISYM1,ISYM12)
          INCX3=NAS1*NISH(ISYM12)
          NFT=NASH(ISYM1)
          NFI=NISH(ISYM1)
          CALL PMLTR1(KOD,IMLTOP,LIST(LLST1),
     &                lg_X,NAS1,NIS1,JXOFF,
     &                WORK(LFTI+IOFFIT(ISYM1)),NFT,NFI,
     &                lg_Y,NAS2,NIS2)
        END IF
C  -----------------------------------------------
      ELSE IF (KOD.EQ.23) THEN
C ICASE1= 5
C ICASE2= 12

C  D&HP One-el
        IF(ISYM1.EQ.1) THEN
          DO ISYMI=1,NSYM
          NI=NISH(ISYMI)
          IF(NI.NE.0) THEN
            ISYMA=ISYMI
            NA=NSSH(ISYMA)
            IF(NA.NE.0) THEN
              ISYMJ=MUL(ISYMI,ISYM2)
              NJ=NISH(ISYMJ)
              IF(NJ.NE.0) THEN
                ISYMB=ISYMJ
                NB=NSSH(ISYMB)
                IF(NB.NE.0) THEN
                  LLST1=LLIST(ISYMI,ISYM2,14)
                  NLST1=NLIST(ISYMI,ISYM2,14)
                  IF(NLST1.NE.0) THEN
                    VAL1(1)=SQRI2
                    VAL1(2)=1.0D0
                    LLST2=LLIST(ISYMA,ISYM2,16)
                    NLST2=NLIST(ISYMA,ISYM2,16)
                    IF(NLST2.NE.0) THEN
                      VAL2(1)=SQRI2
                      VAL2(2)=1.0D0
                      IXIA=1+IOFFIA(ISYMI)
                      INCX1=1
                      INCX2=NI
                      CALL PMLTSCA(KOD,IMLTOP,
     &                             LIST(LLST1),LIST(LLST2),
     &                             X1(IXIA),NI,NA,
     &                             WORK(LFIA+IOFFIA(ISYMB)),NJ,NB,
     &                             lg_Y,NAS2,NIS2)
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END IF
          END DO
        END IF
C ---------------------------
      ELSE IF (KOD.EQ.24) THEN
C ICASE1= 5
C ICASE2= 13

C  D&HM One-el
        IF(ISYM1.EQ.1) THEN
          DO ISYMI=1,NSYM
          NI=NISH(ISYMI)
          IF(NI.NE.0) THEN
            ISYMA=ISYMI
            NA=NSSH(ISYMA)
            IF(NA.NE.0) THEN
              ISYMJ=MUL(ISYMI,ISYM2)
              NJ=NISH(ISYMJ)
              IF(NJ.NE.0) THEN
                ISYMB=ISYMJ
                NB=NSSH(ISYMB)
                IF(NB.NE.0) THEN
                  LLST1=LLIST(ISYMI,ISYM2,15)
                  NLST1=NLIST(ISYMI,ISYM2,15)
                  IF(NLST1.NE.0) THEN
                    VAL1(1)=SQR3*0.5D0
                    VAL1(2)=-VAL1(1)
                    LLST2=LLIST(ISYMA,ISYM2,17)
                    NLST2=NLIST(ISYMA,ISYM2,17)
                    IF(NLST2.NE.0) THEN
                      VAL2(1)=1.0D0
                      VAL2(2)=-1.0D0
                      IXIA=1+IOFFIA(ISYMI)
                      INCX1=1
                      INCX2=NI
                      CALL PMLTSCA(KOD,IMLTOP,
     &                             LIST(LLST1),LIST(LLST2),
     &                             X1(IXIA),NI,NA,
     &                             WORK(LFIA+IOFFIA(ISYMB)),NJ,NB,
     &                             lg_Y,NAS2,NIS2)
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END IF
          END DO
        END IF
C ---------------------------
      ELSE
        WRITE(6,*)' INTERNAL ERROR: SGM reached invalid KOD=',KOD
        Call Abend()
      END IF

      RETURN
      END
