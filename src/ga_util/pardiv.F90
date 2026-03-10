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

function iParDiv(nMax,nMin)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, nProcs
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iParDiv
integer(kind=iwp), intent(in) :: nMax, nMin

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  iParDiv = nMax/nProcs+1+nMin
else
  iParDiv = nMax+nMin
end if
#else
iParDiv = nMax+nMin
#endif

end function iParDiv
