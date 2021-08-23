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

subroutine Get_OrbE(OrbE,nOrbE)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nOrbE
real(kind=wp), intent(out) :: OrbE(nOrbE)
integer(kind=iwp) :: nData
character(len=24) :: Label
logical(kind=iwp) :: Found

Label = 'OrbE'
call qpg_dArray(Label,Found,nData)
if ((.not. Found) .or. (nOrbE == 0)) call SysAbendMsg('get_orbe','Did not find:',Label)
if (nOrbE /= nData) call SysAbendMsg('get_orbe','nOrbE /= nData','')
call Get_dArray(Label,OrbE,nOrbE)

return

end subroutine Get_OrbE
