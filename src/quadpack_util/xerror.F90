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

subroutine xerror(Label,ix,ier,lvl)

use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: ix, ier, lvl

write(u6,*) 'Terminate in xerror!'
write(u6,'(a)') Label
write(u6,'(a,i5)') 'ix=',ix
write(u6,'(a,i5)') 'ier=',ier
write(u6,'(a,i5)') 'lvl=',lvl
call Abend()

return

end subroutine xerror
