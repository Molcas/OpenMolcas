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
      SUBROUTINE TRDNS2D(IVEC,JVEC,DPT2,NDPT2)

#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      IMPLICIT REAL*8 (A-H,O-Z)


#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
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
          CALL RHS_ALLO(NIN,NIS,lg_V1)
          CALL RHS_READ_SR(lg_V1,ICASE,ISYM,IVEC)
          IF(IVEC.EQ.JVEC) THEN
            lg_V2=lg_V1
          ELSE
            CALL RHS_ALLO(NIN,NIS,lg_V2)
            CALL RHS_READ_SR(lg_V2,ICASE,ISYM,JVEC)
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

          CALL RHS_FREE(NIN,NIS,lg_V1)
          IF(IVEC.NE.JVEC) THEN
            CALL RHS_FREE(NIN,NIS,lg_V2)
          END IF
 100    CONTINUE
 101  CONTINUE

      RETURN
      END
