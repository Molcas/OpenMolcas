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
      SUBROUTINE H0SPCT()
      use definitions, only: iwp, wp, u6
      use caspt2_global, only:iPrGlb
      use caspt2_global, only:dnmThr,cntThr,cmpThr
      use caspt2_global, only:LUSBT
      use PrintLevel, only: verbose
#ifdef _MOLCAS_MPP_
      use allgather_wrapper, only : allgather_R, allgather_I
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      use fake_GA, only: GA_Arrays
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      CHARACTER(LEN=80) LINE
      INTEGER, ALLOCATABLE, TARGET:: IDXBUF(:,:)
      REAL*8, ALLOCATABLE, TARGET:: VALBUF(:,:)
#ifdef _MOLCAS_MPP_
      INTEGER, ALLOCATABLE, TARGET:: IDX_H(:,:)
      REAL*8, ALLOCATABLE, TARGET:: VAL_H(:,:)
#endif
      INTEGER, POINTER:: IDX(:,:)=>Null()
      REAL*8, POINTER:: VAL(:,:)=>Null()
      REAL*8, ALLOCATABLE:: BD(:), ID(:)

C Write pertinent warnings and statistics for the energy
C denominators, i.e. the spectrum of (H0(diag)-E0).


      WRITE(u6,*)
      Call CollapseOutput(1,'Denominators, etc.')
      WRITE(u6,'(10A11)')('-----------',i=1,10)
      WRITE(u6,'(A)')' Report on small energy denominators, large'//
     &   ' coefficients, and large energy contributions.'

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(u6,'(A)')
     &   '  The ACTIVE-MIX index denotes linear combinations'//
     &   ' which gives ON expansion functions'
        WRITE(u6,'(A)')'  and makes H0 diagonal within type.'
        WRITE(u6,'(A)')
     &   '  DENOMINATOR: The (H0_ii - E0) value from the above-'//
     &   'mentioned diagonal approximation.'
        WRITE(u6,'(A)')'  RHS VALUE  : Right-Hand Side of CASPT2 Eqs.'
        WRITE(u6,'(A)')
     &   '  COEFFICIENT: Multiplies each of the above ON terms'//
     &   ' in the first-order wave function.'
        WRITE(u6,'(A)')' Thresholds used:'
        WRITE(u6,'(a,f7.4)')'         Denominators:',DNMTHR
        WRITE(u6,'(a,f7.4)')'         Coefficients:',CMPTHR
        WRITE(u6,'(a,f7.4)')' Energy contributions:',CNTTHR
        WRITE(u6,*)
      END IF

      WRITE(u6,'(A)')'CASE  SYMM ACTIVE-MIX  NON-ACTIVE'
     &            //' INDICES          DENOMINATOR'
     &            //'     RHS VALUE       COEFFICIENT'
     &            //'     CONTRIBUTION'

CSVC: initial buffer size, will be reallocated on the fly
      MAXBUF=1024
      CALL mma_allocate(IDXBUF,2,MAXBUF,LABEL='IDXBUF')
      CALL mma_allocate(VALBUF,4,MAXBUF,LABEL='VALBUF')

C Very long loop over symmetry and case:
      DO ICASE=1,13
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF(NIS.EQ.0) GOTO 100
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) GOTO 100
          LINE(1:12)=CASES(ICASE)//'    '
          WRITE(LINE(10:10),'(i1)') ISYM

C Remember: NIN values in BDIAG, but must read NAS for correct
C positioning.
          CALL mma_allocate(BD,NAS,LABEL='BD')
          CALL mma_allocate(ID,NIS,LABEL='ID')
          JD=IDBMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,2,BD,NAS,JD)
          CALL DDAFILE(LUSBT,2,ID,NIS,JD)

          CALL RHS_ALLO(NIN,NIS,lg_RHS)
          CALL RHS_ALLO(NIN,NIS,lg_VEC)
          CALL RHS_READ_SR(lg_RHS,ICASE,ISYM,IRHS)
          CALL RHS_READ_SR(lg_VEC,ICASE,ISYM,IVECX)
          IBUF=0
#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
* Get the superindex ranges of this process's block. If no elements are
* owned by a process, then ilo=0 and ihi=-1 such that the loop further
* down will just be skipped.
            CALL GA_Sync
            myRank = GA_NodeID()
            CALL GA_Distribution (lg_RHS,myRank,IASTA,IAEND,IISTA,IIEND)
            IF (IASTA.NE.0 .AND. IAEND-IASTA+1.NE.NIN) THEN
              WRITE(u6,*) 'RHSOD: mismatch in range of the superindices'
              CALL AbEnd()
            END IF
* if the block is non-empty, loop over its elements
            IF (IASTA.GT.0 .AND. IISTA.GT.0) THEN
              CALL GA_Access (lg_RHS,IASTA,IAEND,IISTA,IIEND,mRHS,LD)
              CALL GA_Access (lg_VEC,IASTA,IAEND,IISTA,IIEND,mVEC,LD)
              IF (LD.NE.NIN) THEN
                WRITE(u6,*) 'RHSOD: assumption NAS=LDW wrong, abort'
                CALL AbEnd()
              END IF
              NA=NAS*(IIEND-IISTA+1)
            END IF
          ELSE
#endif
            IASTA=1
            IAEND=NIN
            IISTA=1
            IIEND=NIS
#ifdef _MOLCAS_MPP_
          END IF
#endif

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
          DO IIS=IISTA,IIEND
            DO IAS=IASTA,IAEND
              DNOM=BD(IAS)+ID(IIS)
#ifdef _MOLCAS_MPP_
              IF (Is_Real_Par()) THEN
                RHS =DBL_MB(mRHS+IAS-1+NIN*(IIS-IISTA))
                COEF=DBL_MB(mVEC+IAS-1+NIN*(IIS-IISTA))
              ELSE
#endif
                RHS =GA_Arrays(lg_RHS)%A(IAS+NIN*(IIS-IISTA))
                COEF=GA_Arrays(lg_VEC)%A(IAS+NIN*(IIS-IISTA))
#ifdef _MOLCAS_MPP_
              END IF
#endif
              ECNT=COEF*RHS
              IF (ABS(DNOM).LT.DNMTHR .OR.
     &            ABS(COEF).GT.CMPTHR .OR.
     &            ABS(ECNT).GT.CNTTHR )
     &        THEN
                IF (IBUF.LT.MAXBUF) THEN
                  IBUF=IBUF+1
                  IDXBUF(1,IBUF)=IAS
                  IDXBUF(2,IBUF)=IIS
                  VALBUF(1,IBUF)=DNOM
                  VALBUF(2,IBUF)=RHS
                  VALBUF(3,IBUF)=COEF
                  VALBUF(4,IBUF)=ECNT
                END IF
              END IF
            END DO
          END DO

#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            IF (IASTA.GT.0 .AND. IISTA.GT.0) THEN
              CALL GA_Release (lg_RHS,IASTA,IAEND,IISTA,IIEND)
              CALL GA_Release (lg_VEC,IASTA,IAEND,IISTA,IIEND)
            END IF
          END IF
#endif
          NBUF=IBUF
#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            CALL GAIGOP_SCAL(NBUF,'+')
            CALL mma_allocate(IDX_H,2,NBUF,LABEL='IDX_H')
            CALL mma_allocate(VAL_H,4,NBUF,LABEL='VAL_H')
            CALL allgather_I(IDXBUF,2*IBUF,IDX_H,2*NBUF)
            CALL allgather_R(VALBUF,4*IBUF,VAL_H,4*NBUF)
            IDX=>IDX_H
            VAL=>VAL_H
          ELSE
#endif
            IDX=>IDXBUF
            VAL=>VALBUF
#ifdef _MOLCAS_MPP_
          END IF
#endif

          DO IBUF=1,NBUF
            IAS  = IDX(1,IBUF)
            IIS  = IDX(2,IBUF)
            DNOM = VAL(1,IBUF)
            RHS  = VAL(2,IBUF)
            COEF = VAL(3,IBUF)
            ECNT = VAL(4,IBUF)
            IF(ICASE.EQ.12.OR.ICASE.EQ.13) THEN
              CALL EXCIND(IAS,IIS,ISYM,ICASE,IP,IQ,IR,IS)
              LINE(13:20)=ORBNAM(IP)
              LINE(21:28)=ORBNAM(IQ)
              LINE(29:36)=ORBNAM(IR)
              LINE(37:44)=ORBNAM(IS)
              LINE(45:46)='  '
            ELSE
              WRITE(LINE(13:22),'(A2,I1,A1,I4.4)')
     &                               'Mu',ISYM,'.',IAS
              CALL NSIND(IIS,ISYM,ICASE,IP,IQ,IR)
              LINE(23:30)=ORBNAM(IP)
              LINE(31:46)='                '
              IF(IQ.GT.0) LINE(31:38)=ORBNAM(IQ)
              IF(IR.GT.0) LINE(39:46)=ORBNAM(IR)
            END IF
            WRITE(u6,'(A,4F16.8)') LINE(1:46),DNOM,RHS,COEF,ECNT
          END DO

#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            CALL mma_deallocate(IDX_H)
            CALL mma_deallocate(VAL_H)
          END IF
#endif

          CALL RHS_FREE(lg_RHS)
          CALL RHS_FREE(lg_VEC)

          CALL mma_deallocate(BD)
          CALL mma_deallocate(ID)

 100      CONTINUE

C End of very long loop over symmetry and case:
        END DO
      END DO

      CALL mma_deallocate(IDXBUF)
      IDX=>Null()
      CALL mma_deallocate(VALBUF)
      VAL=>Null()

      Call CollapseOutput(0,'Denominators, etc.')

      END SUBROUTINE H0SPCT
