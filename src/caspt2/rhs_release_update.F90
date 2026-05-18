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

      SUBROUTINE RHS_RELEASE_UPDATE (lg_W,iLo,iHi,jLo,jHi)
!SVC: this routine releases a local block that was written to back to
! the global array
      use definitions, only: iwp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT None
      integer(kind=iwp), Intent(inout):: lg_W,iLo,iHi,jLo,jHi
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"

      IF (Is_Real_Par()) THEN
        IF (iLo.GT.0 .AND. jLo.GT.0) THEN
          CALL GA_Release_Update (lg_W,iLo,iHi,jLo,jHi)
        END IF
      END IF
#else
#include "macros.fh"
! Avoid unused argument warnings
      unused_var(lg_W)
      unused_var(iLo)
      unused_var(iHi)
      unused_var(jLo)
      unused_var(jHi)
#endif

      END SUBROUTINE RHS_RELEASE_UPDATE
