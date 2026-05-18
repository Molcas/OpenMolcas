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

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHS_INIT()
      use definitions, only: iwp, wp
      use caspt2_global, only: LURHS
      use caspt2_module, only: nSym, nISup, nASup, iOffRHS
      IMPLICIT None
      REAL(kind=wp) DUMMY(1)
      integer(kind=iwp) iDisk, iCase, iSym, NAS, NIS, NW, NRHS, iLo,    &
     &                  iHi, jLo, jHi

!-SVC: loop over symmetry/cases, get local patch of RHS, write, and then
! update the disk address in IOFFRHS
      IDISK=0
      DO ICASE=1,13
        DO ISYM=1,NSYM
          IOFFRHS(ISYM,ICASE)=IDISK

          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          NW=NAS*NIS

          IF (NW.EQ.0) CYCLE

          CALL RHS_DISTRIBUTION(NAS,NIS,iLo,iHi,jLo,jHi)
          NRHS=NAS*(jHi-jLo+1)
          CALL DDAFILE(LURHS(1),0,DUMMY,NRHS,IDISK)

        END DO
      END DO

      END SUBROUTINE RHS_INIT
