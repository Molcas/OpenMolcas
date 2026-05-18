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

subroutine RHS_READ_C(lg_W,iCASE,iSYM,iVEC)

use caspt2_module, only: NASUP, NISUP
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lg_W, iCASE, iSYM, iVEC
integer(kind=iwp) :: NAS, NIS

NAS = NASUP(ISYM,ICASE)
NIS = NISUP(ISYM,ICASE)
call RHS_READ(NAS,NIS,lg_W,ICASE,ISYM,IVEC)

end subroutine RHS_READ_C
