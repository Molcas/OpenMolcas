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

subroutine Get_OrbE_MpProp(ipOrbE,nOrbE)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: ipOrbE, nOrbE
#include "WrkSpc.fh"
character(len=24) :: Label
logical(kind=iwp) :: Found

Label = 'OrbE'
call qpg_dArray(Label,Found,nOrbE)
if ((.not. Found) .or. (nOrbE == 0)) then
  call SysAbendMsg('get_orbe','Did not find:',Label)
end if
call GetMem('OrbE','Allo','Real',ipOrbE,2*nOrbE)
call Get_dArray(Label,Work(ipOrbE),nOrbE)

return

end subroutine Get_OrbE_MpProp
