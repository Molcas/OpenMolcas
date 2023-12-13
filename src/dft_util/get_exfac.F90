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

function Get_ExFac(KSDFT)
!***********************************************************************
! Return the factor which determines how much "exact exchange" that    *
! should be included.                                                  *
!***********************************************************************

use Functionals, only: Get_Func_ExFac
use Constants, only: Zero, One
use Definitions, only: wp

implicit none
real(kind=wp) :: Get_ExFac
character(len=*), intent(in) :: KSDFT
character(len=80) :: cTmp

!                                                                      *
!***********************************************************************
!                                                                      *
! Write functional to run file.

if (KSDFT /= 'Overlap') then
  cTmp = KSDFT
  call Put_cArray('DFT functional',cTmp,80)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if ((KSDFT(1:2) == 'T:') .or. (KSDFT(1:3) == 'FT:')) then
  Get_ExFac = Zero
  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! We bring in only cases where it is different from zero.
select case (KSDFT)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  case ('CASDFT')
    Get_ExFac = One
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  case ('SCF')
    Get_ExFac = One
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  case ('CS')
    Get_ExFac = One
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  case default
    Get_ExFac = Get_Func_ExFac(KSDFT)
!                                                                      *
!***********************************************************************
!                                                                      *
end select
!                                                                      *
!***********************************************************************
!                                                                      *

end function Get_ExFac
