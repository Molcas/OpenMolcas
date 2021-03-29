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

subroutine BasisConsistency(LuWr,iErr)

use ZMatConv_Mod, only: BasAva, BasReq
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: LuWr
integer(kind=iwp), intent(out) :: iErr
integer(kind=iwp) :: i

iErr = 0
do i=1,size(BasReq)
  if (BasReq(i) .and. (.not. BasAva(i))) then
    iErr = 1
    write(LuWr,*) ' [BasisConsistency]: Atom NA=',i,' requires BS'
    exit
  end if
end do

return

end subroutine BasisConsistency
