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
! Written by Steven Vancoillie, May 2011
! A set of subroutines that can handle RHS arrays in either a serial or
! parallel environment, depending on the situation.
!***********************************************************************
! --> when running serially, the RHS arrays are stored on LUSOLV and are
! loaded into the WORK array when needed.
! --> when running in parallel, the RHS arrays are stored on disk as
! disk resident arrays (DRAs) with filename RHS_XX_XX_XX, where XX is a
! number referring to the case, symmetry, and RHS vector respectively,
! and are loaded onto a global array when needed.
!***********************************************************************

      SUBROUTINE RHS_FPRINT(CTYPE,IVEC)
      use definitions, only: iwp, wp, u6
      use constants, only: Zero
      use caspt2_module, only: NSYM, NASUP, NINDEP, NISUP
      IMPLICIT None

      integer(kind=iwp), Intent(in):: IVEC
      CHARACTER(LEN=*), intent(in):: CTYPE

      real(kind=wp) :: FP(8)
      integer(kind=iwp) NROW, ICASE, ISYM, NAS, NIN, NIS, lg_W
      real(kind=wp), external:: RHS_DDOT

!-SVC: print out DNRM2 of the all RHS components
      NROW=0 ! dummy initialize
      DO ICASE=1,13
        DO ISYM=1,NSYM

          NAS=NASUP(ISYM,ICASE)
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)

          IF (CTYPE.EQ.'C') THEN
            NROW=NAS
          ELSE IF (CTYPE.EQ.'SR') THEN
            NROW=NIN
          ELSE
            WRITE(u6,'(1X,A)') 'RHS_FPRINT: invalid type: '//CTYPE
            CALL ABEND()
          END IF

          IF (NAS.NE.0 .AND. NIN.NE.0 .AND. NIS.NE.0) THEN
            CALL RHS_ALLO(NROW,NIS,lg_W)
            CALL RHS_READ(NROW,NIS,lg_W,iCASE,iSYM,iVEC)
            FP(ISYM)=SQRT(RHS_DDOT(NROW,NIS,lg_W,lg_W))
            CALL RHS_FREE(lg_W)
          ELSE
            FP(ISYM)=Zero
          END IF
        END DO
        WRITE(u6,'(1X,I2,1X,8F21.14)') ICASE, FP(1:NSYM)
      END DO

      END SUBROUTINE RHS_FPRINT
