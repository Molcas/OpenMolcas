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

subroutine Chk_LblCnt(Lbl,nList)

use Center_Info, only: dc
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: Lbl
integer(kind=iwp), intent(in) :: nList
integer(kind=iwp) :: iList
character(len=72) :: Warning

do iList=1,nList
  if (Lbl == dc(iList)%LblCnt) then
    write(Warning,'(A,A)') 'ChkLbl: Duplicate label; Lbl=',Lbl
    call WarningMessage(2,Warning)
    call Quit_OnUserError()
  end if
end do

return

end subroutine Chk_LblCnt
