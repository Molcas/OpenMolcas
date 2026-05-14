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

subroutine Mk_Indkl(Indkl_OnOff,Indkl,nkl)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nkl, Indkl_OnOff(nkl)
integer(kind=iwp), intent(out) :: Indkl(nkl)
integer(kind=iwp) :: ikl, jkl

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call iVcPrt('Mk_Indkl: Indkl_OnOff',' ',Indkl_OnOff,nkl)
#endif
ikl = 0
do jkl=1,nkl
  if (Indkl_OnOff(jkl) == 1) then
    ikl = ikl+1
    Indkl(jkl) = ikl
  else
    Indkl(jkl) = 0
  end if
end do
#ifdef _DEBUGPRINT_
call iVcPrt('Mk_Indkl: Indkl',' ',Indkl,nkl)
#endif

return

end subroutine Mk_Indkl
