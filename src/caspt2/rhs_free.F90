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

subroutine RHS_FREE(lg_W)
!SVC: this routine writes the RHS array to disk

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use fake_GA, only: Deallocate_GA_Array
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: lg_w
#ifdef _MOLCAS_MPP_
logical(kind=iwp) :: bStat
#include "global.fh"
#include "mafdecls.fh"

if (Is_Real_Par()) then
  !SVC: Destroy the global array
  bStat = GA_Destroy(lg_W)
# include "macros.fh"
  unused_var(bStat)
else
#endif
  call Deallocate_GA_Array(lg_W)
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine RHS_FREE
