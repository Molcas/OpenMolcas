!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2011, Steven Vancoillie                                *
!***********************************************************************
!***********************************************************************
! New VECTOR UTILITIES, written by Steven Vancoillie, May 2011
! A set of subroutines that can transform RHS arrays using the parallel
! aware subroutines
!***********************************************************************

      SUBROUTINE PTRTOSR (ITYPE,IVEC,JVEC)
      use Constants, only: Zero
      use caspt2_module, only: nSym, nInDep, nASup, nISup, nCases,      &
     &                         CPUVec, TIOVec
      use definitions, only: iwp, wp
      IMPLICIT None

      integer(kind=iwp), intent(In):: iTYPE, iVec, jVec

      integer(kind=iwp) :: iCase, iSym, nIn, nAs, nIs
      integer(kind=iwp) :: lg_v1, lg_v2
      real(kind=wp) CPU1, CPU0, TIO1, TIO0, CPU, TIO


! Transform RHS vectors from SR format to C format.
! ITYPE=0 uses only T matrix, ITYPE=1 uses S*T matrix

      CALL TIMING(CPU0,CPU,TIO0,TIO)

      DO ICASE=1,NCASES
        IF(IVEC.EQ.JVEC .AND. ICASE.EQ.12) Cycle
        IF(IVEC.EQ.JVEC .AND. ICASE.EQ.13) Cycle
        DO ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) Cycle
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF(NIS.EQ.0) Cycle
          CALL RHS_ALLO (NIN,NIS,lg_V1)
          IF(ICASE.NE.12 .AND. ICASE.NE.13) THEN
            IF(NAS.GT.0) THEN
              CALL RHS_ALLO (NAS,NIS,lg_V2)
              CALL RHS_READ (NAS,NIS,lg_V2,ICASE,ISYM,IVEC)
              CALL RHS_SR2C (ITYPE,1,NAS,NIS,NIN,                       &
     &                       lg_V1,lg_V2,ICASE,ISYM)
              CALL RHS_FREE (lg_V2)
            ELSE
              CALL RHS_SCAL (NIN,NIS,lg_V1,Zero)
            END IF
          ELSE
            CALL RHS_READ (NIN,NIS,lg_V1,ICASE,ISYM,IVEC)
          END IF

          CALL RHS_SAVE (NIN,NIS,lg_V1,ICASE,ISYM,JVEC)

          CALL RHS_FREE (lg_V1)

        End Do
      End Do

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUVEC=CPUVEC+(CPU1-CPU0)
      TIOVEC=TIOVEC+(TIO1-TIO0)

      END SUBROUTINE PTRTOSR
