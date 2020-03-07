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
      SUBROUTINE TRDNS1(IVEC,DPT1)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "sigma.fh"
      DIMENSION DPT1(*)
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      CALL QENTER('TRDNS1')
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
      CALL GETMEM('WTI','ALLO','REAL',LWTI,NWTI)
      CALL DCOPY_(NWTI,[0.0D0],0,WORK(LWTI),1)
      ICASE=1
      IWOFF=0
      DO 100 ISYM=1,NSYM
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 100
        NIS=NISUP(ISYM,ICASE)
        NAS=NASUP(ISYM,ICASE)
        NVEC=NIS*NAS
        IF(NVEC.EQ.0) GOTO 100
        CALL RHS_ALLO(NAS,NIS,LVEC)
        CALL RHS_READ_C(LVEC,ICASE,ISYM,IVEC)
        !CALL GETMEM('VEC','ALLO','REAL',LVEC,NVEC)
        !CALL RDBLKC(ISYM,ICASE,IVEC,WORK(LVEC))
        FACT=1.0D00/(DBLE(MAX(1,NACTEL)))
#ifdef _MOLCAS_MPP_
        IF (IS_REAL_PAR()) THEN
          IF (KING()) THEN
            CALL GETMEM('TMP','ALLO','REAL',LTMP,NVEC)
            CALL RHS_GET (NAS,NIS,LVEC,WORK(LTMP))
            CALL SPEC1A(IMLTOP,FACT,ISYM,
     &                  WORK(LTMP),WORK(LWTI+IWOFF))
            CALL GETMEM('TMP','FREE','REAL',LTMP,NVEC)
          END IF
        ELSE
          CALL SPEC1A(IMLTOP,FACT,ISYM,WORK(LVEC),
     &                                 WORK(LWTI+IWOFF))
        END IF
#else
        CALL SPEC1A(IMLTOP,FACT,ISYM,WORK(LVEC),
     &                               WORK(LWTI+IWOFF))
#endif
        !CALL GETMEM('VEC','FREE','REAL',LVEC,NVEC)
        CALL RHS_FREE(NAS,NIS,LVEC)
        IWOFF=IWOFF+NASH(ISYM)*NISH(ISYM)
 100  CONTINUE
 110  CONTINUE

      IF(NWAT.EQ.0) GOTO 210
      CALL GETMEM('WAT','ALLO','REAL',LWAT,NWAT)
      CALL DCOPY_(NWAT,[0.0D0],0,WORK(LWAT),1)
      ICASE=4
      IWOFF=0
      DO 200 ISYM=1,NSYM
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 200
        NIS=NISUP(ISYM,ICASE)
        NAS=NASUP(ISYM,ICASE)
        NVEC=NIS*NAS
        IF(NVEC.EQ.0) GOTO 200
        !CALL GETMEM('VEC','ALLO','REAL',LVEC,NVEC)
        !CALL RDBLKC(ISYM,ICASE,IVEC,WORK(LVEC))
        CALL RHS_ALLO(NAS,NIS,LVEC)
        CALL RHS_READ_C(LVEC,ICASE,ISYM,IVEC)
        FACT=1.0D00/(DBLE(MAX(1,NACTEL)))
#ifdef _MOLCAS_MPP_
        IF (IS_REAL_PAR()) THEN
          IF (KING()) THEN
            CALL GETMEM('TMP','ALLO','REAL',LTMP,NVEC)
            CALL RHS_GET (NAS,NIS,LVEC,WORK(LTMP))
            CALL SPEC1C(IMLTOP,FACT,ISYM,
     &                  WORK(LTMP),WORK(LWAT+IWOFF))
            CALL GETMEM('TMP','FREE','REAL',LTMP,NVEC)
          END IF
        ELSE
          CALL SPEC1C(IMLTOP,FACT,ISYM,WORK(LVEC),
     &                                 WORK(LWAT+IWOFF))
        END IF
#else
        CALL SPEC1C(IMLTOP,FACT,ISYM,WORK(LVEC),
     &                               WORK(LWAT+IWOFF))
#endif
        !CALL GETMEM('VEC','FREE','REAL',LVEC,NVEC)
        CALL RHS_FREE(NAS,NIS,LVEC)
        IWOFF=IWOFF+NSSH(ISYM)*NASH(ISYM)
 200  CONTINUE
 210  CONTINUE

      IF(NWAI.EQ.0) GOTO 300
      CALL GETMEM('WAI','ALLO','REAL',LWAI,NWAI)
      CALL DCOPY_(NWAI,[0.0D0],0,WORK(LWAI),1)
      ICASE=5
      ISYM=1
      IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 300
      NIS=NISUP(ISYM,ICASE)
      NAS=NASUP(ISYM,ICASE)
      NVEC=NIS*NAS
      IF(NVEC.EQ.0) GOTO 300
      !CALL GETMEM('VEC','ALLO','REAL',LVEC,NVEC)
      !CALL RDBLKC(ISYM,ICASE,IVEC,WORK(LVEC))
      CALL RHS_ALLO(NAS,NIS,LVEC)
      CALL RHS_READ_C(LVEC,ICASE,ISYM,IVEC)
      FACT=1.0D00/(DBLE(MAX(1,NACTEL)))
#ifdef _MOLCAS_MPP_
      IF (IS_REAL_PAR()) THEN
        IF (KING()) THEN
          CALL GETMEM('TMP','ALLO','REAL',LTMP,NVEC)
          CALL RHS_GET (NAS,NIS,LVEC,WORK(LTMP))
          CALL SPEC1D(IMLTOP,FACT,WORK(LTMP),WORK(LWAI))
          CALL GETMEM('TMP','FREE','REAL',LTMP,NVEC)
        END IF
      ELSE
        CALL SPEC1D(IMLTOP,FACT,WORK(LVEC),WORK(LWAI))
      END IF
#else
      CALL SPEC1D(IMLTOP,FACT,WORK(LVEC),WORK(LWAI))
#endif
      !CALL GETMEM('VEC','FREE','REAL',LVEC,NVEC)
      CALL RHS_FREE(NAS,NIS,LVEC)
 300  CONTINUE

C Transform vectors back to eigenbasis of H0(diag).
      CALL PTRTOSR(0,IVEC,IVEC)

      IF(NWTI.GT.0) CALL GADSUM(WORK(LWTI),NWTI)
      IF(NWAI.GT.0) CALL GADSUM(WORK(LWAI),NWAI)
      IF(NWAT.GT.0) CALL GADSUM(WORK(LWAT),NWAT)
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
            DPT1(ID)=DPT1(ID)+WORK(LWTI-1+IW)
          END DO
        END DO
        DO IT=1,NA
          ITTOT=IT+NI
          ID=IDOFF+NI+NA+NO*(ITTOT-1)
          IW=IWAT+IT-NA
          DO IA=1,NS
            ID=ID+1
            IW=IW+NA
            DPT1(ID)=DPT1(ID)+WORK(LWAT-1+IW)
          END DO
        END DO
        DO II=1,NI
          ID=IDOFF+NI+NA+NO*(II-1)
          IW=IWAI+II-NI
          DO IA=1,NS
            ID=ID+1
            IW=IW+NI
            DPT1(ID)=DPT1(ID)+WORK(LWAI-1+IW)
          END DO
        END DO
        IWTI=IWTI+NA*NI
        IWAI=IWAI+NS*NI
        IWAT=IWAT+NS*NA
        IDOFF=IDOFF+NO**2
      END DO

      IF(NWTI.GT.0)CALL GETMEM('WTI','FREE','REAL',LWTI,NWTI)
      IF(NWAI.GT.0)CALL GETMEM('WAI','FREE','REAL',LWAI,NWAI)
      IF(NWAT.GT.0)CALL GETMEM('WAT','FREE','REAL',LWAT,NWAT)

      CALL QEXIT('TRDNS1')
      RETURN
      END
