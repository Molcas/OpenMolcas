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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine loopstr0_cvb(iocc,indx,nel,norb)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nel, norb
integer(kind=iwp), intent(out) :: iocc(nel), indx
integer(kind=iwp) :: iel

if (nel > norb) then
  write(u6,*) ' Error in loopstr0, nel > norb :',nel,norb
  call abend_cvb()
end if
indx = 1
do iel=1,nel
  iocc(iel) = iel
end do

return

end subroutine loopstr0_cvb
