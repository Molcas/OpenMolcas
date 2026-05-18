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

      SUBROUTINE PSCAVEC (FACT,IVEC,JVEC)
      use constants, only: Zero, One
      use caspt2_global, ONLY: iPrGlb
      USE PrintLevel, ONLY: USUAL
      use caspt2_module, only: CPUSCA, nCases, nSym, TIOSCA, nInDep,    &
     &                         niSup
      use definitions, only: iwp, wp, u6

      IMPLICIT NONE

      real(kind=wp), intent(in) ::  FACT
      integer(kind=iwp), intent(in) :: IVEC,JVEC

      real(kind=wp) :: CPU0,CPU1,CPU,TIO0,TIO1,TIO
      integer(kind=iwp) :: ICASE,ISYM,NIN,NIS
      integer(kind=iwp) :: lg_V

      REAL(kind=wp) ::  SIGMA2
      real(kind=wp), EXTERNAL :: RHS_DDOT

! Scale vector nr IVEC with scale factor FACT and put the result in
! vector nr JVEC: |JVEC> <- FACT * |IVEC>
      CALL TIMING(CPU0,CPU,TIO0,TIO)

      IF(FACT.EQ.One.AND.IVEC.EQ.JVEC) RETURN
      SIGMA2=Zero
      DO ICASE=1,NCASES
        DO ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF (NIN*NIS.NE.0) THEN
            CALL RHS_ALLO (NIN,NIS,lg_V)
            CALL RHS_READ (NIN,NIS,lg_V,ICASE,ISYM,IVEC)
            CALL RHS_SCAL (NIN,NIS,lg_V,FACT)
            IF(FACT.EQ.-One)SIGMA2=SIGMA2                               &
     &                               +RHS_DDOT(NIN,NIS,lg_V,lg_V)
            CALL RHS_SAVE (NIN,NIS,lg_V,ICASE,ISYM,JVEC)
            CALL RHS_FREE (lg_V)
          END IF
        END DO
      END DO
      IF ((IPRGLB.GE.USUAL).AND.(FACT.EQ.-One)) THEN ! it is at ITER=0
         WRITE(u6,*)
         WRITE(u6,'(1x,a,f18.10)') 'Variance of |WF0>: ',SIGMA2
      END IF

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUSCA=CPUSCA+(CPU1-CPU0)
      TIOSCA=TIOSCA+(TIO1-TIO0)

      END SUBROUTINE PSCAVEC
