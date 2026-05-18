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

subroutine RHS_SCAL(NAS,NIS,lg_W,FACT)
!SVC: this routine multiplies the RHS array with FACT

use definitions, only: iwp, wp
use constants, only: Zero, One
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use fake_GA, only: GA_Arrays

implicit none
integer(kind=iwp), intent(in) :: NAS, NIS, lg_W
real(kind=wp), intent(in) :: FACT
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  if (FACT == Zero) then
    !call GA_Fill (lg_W,Zero)
    call GA_Zero(lg_W)
  else if (FACT /= One) then
    call GA_Scale(lg_W,FACT)
  end if
else
#endif
  if (FACT == Zero) then
    call DCOPY_(NAS*NIS,[Zero],0,GA_Arrays(lg_W)%A,1)
  else if (FACT /= One) then
    call DSCAL_(NAS*NIS,FACT,GA_Arrays(lg_W)%A,1)
  end if
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine RHS_SCAL
