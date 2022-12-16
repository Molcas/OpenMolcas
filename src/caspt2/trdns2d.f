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
      SUBROUTINE TRDNS2D(IVEC,JVEC,DPT2,NDPT2,SCAL)

      use caspt2_global, only: imag_shift
      use caspt2_gradient, only: do_grad
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      IMPLICIT REAL*8 (A-H,O-Z)


#include "rasdim.fh"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "sigma.fh"

      DIMENSION DPT2(NDPT2)

C Add to the diagonal blocks of transition density matrix,
C    DPT2(p,q) = Add <IVEC| E(p,q) |JVEC>,
C i.e. inactive/inactive, active/active, and virt/virt
C submatrices.IVEC, JVEC stands for the 1st-order perturbed
C CASPT2 wave functions in vectors nr IVEC, JVEC on LUSOLV.

C Inact/Inact and Virt/Virt blocks:
      DO 101 ICASE=1,13
        DO 100 ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) GOTO 100
          NIS=NISUP(ISYM,ICASE)
          NVEC=NIN*NIS
          IF(NVEC.EQ.0) GOTO 100
          !! lg_V1: T+lambda
          !! lg_V2: T
          !! IVEC = iVecX
          !! JVEC = iVecR
          CALL RHS_ALLO(NIN,NIS,lg_V1)
          CALL RHS_READ_SR(lg_V1,ICASE,ISYM,IVEC)
          IF(IVEC.EQ.JVEC) THEN
            lg_V2=lg_V1
          ELSE
            CALL RHS_ALLO(NIN,NIS,lg_V2)
            CALL RHS_READ_SR(lg_V2,ICASE,ISYM,JVEC)
            If (do_grad) Then
              If (Scal.ne.1.0D+00)
     *          Call DScal_(NIN*NIS,Scal,Work(lg_V1),1)
            Call DaXpY_(nIN*nIS,1.0D+00,Work(lg_V2),1,Work(lg_V1),1)
            CALL RHS_READ_SR(lg_V2,ICASE,ISYM,IVEC)
            End If
          END IF

CSVC: DIADNS can currently not handle pieces of RHS, so pass the
C full array in case we are running in parallel
#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            IF (KING()) THEN
              ! copy global array to local buffer
              CALL GETMEM('VEC1','ALLO','REAL',LVEC1,NVEC)
              CALL GA_GET(lg_V1,1,NIN,1,NIS,WORK(LVEC1),NIN)
              IF(IVEC.EQ.JVEC) THEN
                LVEC2=LVEC1
              ELSE
                CALL GETMEM('VEC2','ALLO','REAL',LVEC2,NVEC)
                CALL GA_GET(lg_V2,1,NIN,1,NIS,WORK(LVEC2),NIN)
              END IF

              CALL DIADNS(ISYM,ICASE,WORK(LVEC1),WORK(LVEC2),
     &                    DPT2,iWORK(LLISTS))

              ! free local buffer
              CALL GETMEM('VEC1','FREE','REAL',LVEC1,NVEC)
              IF(IVEC.NE.JVEC) THEN
                CALL GETMEM('VEC2','FREE','REAL',LVEC2,NVEC)
              END IF
            END IF
            CALL GASYNC
          ELSE
            CALL DIADNS(ISYM,ICASE,WORK(lg_V1),WORK(lg_V2),
     &                  DPT2,iWORK(LLISTS))
          END IF
#else
          CALL DIADNS(ISYM,ICASE,WORK(lg_V1),WORK(lg_V2),
     &                DPT2,iWORK(LLISTS))
#endif
          If (do_grad .and. imag_shift .ne. 0.0d0) Then
            nAS = nASUP(iSym,iCase)
            Call GETMEM('LBD','ALLO','REAL',LBD,nAS)
            Call GETMEM('LID','ALLO','REAL',LID,nIS)
            iD = iDBMat(iSym,iCase)
            Call dDaFile(LUSBT,2,Work(LBD),nAS,iD)
            Call dDaFile(LUSBT,2,Work(LID),nIS,iD)
C
            CALL RHS_READ_SR(lg_V1,ICASE,ISYM,IVEC)
            Call CASPT2_ResD(2,nIN,nIS,lg_V1,Work(LBD),Work(LID))
            CALL RHS_READ_SR(lg_V2,ICASE,ISYM,JVEC)
            Call CASPT2_ResD(2,nIN,nIS,lg_V2,Work(LBD),Work(LID))
C
            Call DScal_(NDPT2,-1.0D+00,DPT2,1)
#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              IF (KING()) THEN
                ! copy global array to local buffer
                CALL GETMEM('VEC1','ALLO','REAL',LVEC1,NVEC)
                CALL GA_GET(lg_V1,1,NIN,1,NIS,WORK(LVEC1),NIN)
                IF(IVEC.EQ.JVEC) THEN
                  LVEC2=LVEC1
                ELSE
                  CALL GETMEM('VEC2','ALLO','REAL',LVEC2,NVEC)
                  CALL GA_GET(lg_V2,1,NIN,1,NIS,WORK(LVEC2),NIN)
                END IF

                CALL DIADNS(ISYM,ICASE,WORK(LVEC1),WORK(LVEC2),
     &                      DPT2,iWORK(LLISTS))

                ! free local buffer
                CALL GETMEM('VEC1','FREE','REAL',LVEC1,NVEC)
                IF(IVEC.NE.JVEC) THEN
                  CALL GETMEM('VEC2','FREE','REAL',LVEC2,NVEC)
                END IF
              END IF
              CALL GASYNC
            ELSE
              CALL DIADNS(ISYM,ICASE,WORK(lg_V1),WORK(lg_V2),
     &                    DPT2,iWORK(LLISTS))
            END IF
#else
            CALL DIADNS(ISYM,ICASE,WORK(lg_V1),WORK(lg_V2),
     &                  DPT2,iWORK(LLISTS))
#endif
            Call DScal_(NDPT2,-1.0D+00,DPT2,1)
            Call GETMEM('LBD','FREE','REAL',LBD,nAS)
            Call GETMEM('LID','FREE','REAL',LID,nIS)
          End IF

          CALL RHS_FREE(NIN,NIS,lg_V1)
          IF(IVEC.NE.JVEC) THEN
            CALL RHS_FREE(NIN,NIS,lg_V2)
          END IF
 100    CONTINUE
 101  CONTINUE

      RETURN
      END
