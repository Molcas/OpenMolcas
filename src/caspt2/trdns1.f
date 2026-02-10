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
      SUBROUTINE TRDNS1(IVEC,DPT1,NDPT1)
      use definitions, only: iwp, wp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      use stdalloc, only: mma_allocate, mma_deallocate
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: nActel, nSym, nIsh, nAsh, nSsh, nInDep,
     &                         nISup, nASup, nAsh, nOrb
#define RHS_ X_RHS_
      IMPLICIT none

      integer(kind=iwp), intent(in):: IVEC, NDPT1
      real(kind=wp), intent(inout):: DPT1(NDPT1)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      real(kind=wp), ALLOCATABLE:: WTI(:), WAT(:), WAI(:)
#ifdef _MOLCAS_MPP_
      real(kind=wp), ALLOCATABLE:: TMP(:)
#endif
      real(kind=wp) FACT
      integer(kind=iwp) IA, ICASE, ID, IDOFF, II, IMLTOP, ISYM, IT,
     &                  ITTOT, IW, IWAI, IWAT, IWOFF, IWTI, lVec, NA,
     &                  NAS, NI, NIS, NO, NS, NVEC, NWAI, NWAT, NWTI

C Add to the transition density matrix DPT1,
C    DPT1(p,q) = Add <IVEC| E(p,q) |0>.
C where IVEC stands for the 1st-order perturbed CASPT2
C wave function stored as vector nr IVEC on LUSOLV.
C DPT1 is stored as symmetry-blocked array of square matrices.
C Each square matrix is actually lower block triangle, but is
C stored in full, including zero elements.

C Only cases A, C and D(Symm 1) contributes.

C Transform to standard representation, covariant form.
      CALL PTRTOC(1,IVEC,IVEC)

      NWTI=0
      NWAI=0
      NWAT=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NS=NSSH(ISYM)
        NWTI=NWTI+NA*NI
        NWAI=NWAI+NS*NI
        NWAT=NWAT+NS*NA
      END DO

      IMLTOP=1
      IF(NWTI.EQ.0) GOTO 110
      CALL mma_allocate(WTI,NWTI,LABEL='WTI')
      WTI(:)=0.0D0
      ICASE=1
      IWOFF=1
      DO 100 ISYM=1,NSYM
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 100
        NIS=NISUP(ISYM,ICASE)
        NAS=NASUP(ISYM,ICASE)
        NVEC=NIS*NAS
        IF(NVEC.EQ.0) GOTO 100
        CALL RHS_ALLO(NAS,NIS,LVEC)
        CALL RHS_READ_C(LVEC,ICASE,ISYM,IVEC)
        FACT=1.0D00/(DBLE(MAX(1,NACTEL)))
#ifdef _MOLCAS_MPP_
        IF (IS_REAL_PAR()) THEN
          IF (KING()) THEN
            CALL mma_allocate(TMP,NVEC,Label='TMP')
            CALL RHS_GET (NAS,NIS,LVEC,TMP)
            CALL SPEC1A(IMLTOP,FACT,ISYM,TMP,WTI(IWOFF))
            CALL mma_deallocate(TMP)
          END IF
        ELSE
#endif
          CALL SPEC1A(IMLTOP,FACT,ISYM,GA_Arrays(LVEC)%A,
     &                                 WTI(IWOFF))
#ifdef _MOLCAS_MPP_
        END IF
#endif
        CALL RHS_FREE(LVEC)
        IWOFF=IWOFF+NASH(ISYM)*NISH(ISYM)
 100  CONTINUE
 110  CONTINUE

      IF(NWAT.EQ.0) GOTO 210
      CALL mma_allocate(WAT,NWAT,Label='WAT')
      WAT(:)=0.0D0
      ICASE=4
      IWOFF=1
      DO 200 ISYM=1,NSYM
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 200
        NIS=NISUP(ISYM,ICASE)
        NAS=NASUP(ISYM,ICASE)
        NVEC=NIS*NAS
        IF(NVEC.EQ.0) GOTO 200
        IF(NSSH(ISYM)*NASH(ISYM).EQ.0) GOTO 200
        CALL RHS_ALLO(NAS,NIS,LVEC)
        CALL RHS_READ_C(LVEC,ICASE,ISYM,IVEC)
        FACT=1.0D00/(DBLE(MAX(1,NACTEL)))
#ifdef _MOLCAS_MPP_
        IF (IS_REAL_PAR()) THEN
          IF (KING()) THEN
            CALL mma_allocate(TMP,NVEC,LABEL='TMP')
            CALL RHS_GET (NAS,NIS,LVEC,TMP)
            CALL SPEC1C(IMLTOP,FACT,ISYM,TMP,WAT(IWOFF))
            CALL mma_deallocate(TMP)
          END IF
        ELSE
#endif
          CALL SPEC1C(IMLTOP,FACT,ISYM,GA_Arrays(LVEC)%A,WAT(IWOFF))
#ifdef _MOLCAS_MPP_
        END IF
#endif
        CALL RHS_FREE(LVEC)
        IWOFF=IWOFF+NSSH(ISYM)*NASH(ISYM)
 200  CONTINUE
 210  CONTINUE

      IF(NWAI.EQ.0) GOTO 300
      CALL mma_allocate(WAI,NWAI,Label='WAI')
      WAI(:)=0.0D0
      ICASE=5
      ISYM=1
      IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 300
      NIS=NISUP(ISYM,ICASE)
      NAS=NASUP(ISYM,ICASE)
      NVEC=NIS*NAS
      IF(NVEC.EQ.0) GOTO 300
      CALL RHS_ALLO(NAS,NIS,LVEC)
      CALL RHS_READ_C(LVEC,ICASE,ISYM,IVEC)
      FACT=1.0D00/(DBLE(MAX(1,NACTEL)))
#ifdef _MOLCAS_MPP_
      IF (IS_REAL_PAR()) THEN
        IF (KING()) THEN
          CALL mma_allocate(TMP,NVEC,LABEL='TMP')
          CALL RHS_GET (NAS,NIS,LVEC,TMP)
          CALL SPEC1D(IMLTOP,FACT,TMP,WAI)
          CALL mma_deallocate(TMP)
        END IF
      ELSE
#endif
        CALL SPEC1D(IMLTOP,FACT,GA_Arrays(LVEC)%A,WAI)
#ifdef _MOLCAS_MPP_
      END IF
#endif
      CALL RHS_FREE(LVEC)
 300  CONTINUE

C Transform vectors back to eigenbasis of H0(diag).
      CALL PTRTOSR(0,IVEC,IVEC)

      IF(NWTI.GT.0) CALL GADSUM(WTI,NWTI)
      IF(NWAI.GT.0) CALL GADSUM(WAI,NWAI)
      IF(NWAT.GT.0) CALL GADSUM(WAT,NWAT)
C Put transition density elements in temporaries W into
C proper positions, as subdiagonal matrices in DPT1:
      IDOFF=0
      IWTI=0
      IWAI=0
      IWAT=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NS=NSSH(ISYM)
        NO=NORB(ISYM)
        DO IT=1,NA
          ITTOT=IT+NI
          DO II=1,NI
            ID=IDOFF+ITTOT+NO*(II-1)
            IW=IWTI+IT+NA*(II-1)
            DPT1(ID)=DPT1(ID)+WTI(IW)
          END DO
        END DO
        DO IT=1,NA
          ITTOT=IT+NI
          ID=IDOFF+NI+NA+NO*(ITTOT-1)
          IW=IWAT+IT-NA
          DO IA=1,NS
            ID=ID+1
            IW=IW+NA
            DPT1(ID)=DPT1(ID)+WAT(IW)
          END DO
        END DO
        DO II=1,NI
          ID=IDOFF+NI+NA+NO*(II-1)
          IW=IWAI+II-NI
          DO IA=1,NS
            ID=ID+1
            IW=IW+NI
            DPT1(ID)=DPT1(ID)+WAI(IW)
          END DO
        END DO
        IWTI=IWTI+NA*NI
        IWAI=IWAI+NS*NI
        IWAT=IWAT+NS*NA
        IDOFF=IDOFF+NO**2
      END DO

      IF(NWTI.GT.0)CALL mma_deallocate(WTI)
      IF(NWAI.GT.0)CALL mma_deallocate(WAI)
      IF(NWAT.GT.0)CALL mma_deallocate(WAT)

      END SUBROUTINE TRDNS1
