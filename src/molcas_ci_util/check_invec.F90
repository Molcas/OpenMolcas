!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

! PAM2009 Who wrote this? What purpose?
subroutine Check_InVec(InVec)

#ifdef _MOLCAS_MPP_
use Para_Info, only: nProcs, Is_Real_Par
use Definitions, only: u6
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: InVec
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: InVec_Tot, mProcs

if (.not. Is_Real_Par()) return
InVec_Tot = InVec
call GAIGOP_SCAL(InVec_Tot,'+')
if (InVec_Tot /= 0) then
  if (InVec == 0) then
    mProcs = -1
  else
    mProcs = InVec_Tot/InVec
  end if
  if (mProcs /= nProcs) then
    write(u6,*) 'Check_InVec: different orbital options on different nodes'
    write(u6,*) 'Sets InVec to 0'
    InVec = 0
  end if
end if

#else

#include "macros.fh"
unused_var(InVec)

#endif

end subroutine Check_InVec
