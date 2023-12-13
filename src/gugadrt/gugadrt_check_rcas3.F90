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

subroutine gugadrt_check_rcas3(jk,ind,inb,ndj,locu)

use gugadrt_global, only: ja, jb, max_node
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: jk, ind(8,max_node), ndj, locu(8,ndj)
integer(kind=iwp), intent(out) :: inb
integer(kind=iwp) :: i, iex, iexcit(ndj), lsym(8), m, nsumel

inb = 0
nsumel = 0
do i=1,8
  lsym(i) = ind(i,jk)
  nsumel = nsumel+lsym(i)
end do
do i=1,ndj
  iexcit(i) = 0
  do m=1,8
    iex = lsym(m)-locu(m,i)
    if (iex > 0) then
      iexcit(i) = iexcit(i)+iex
    end if
  end do
end do
inb = minval(iexcit)+ja(jk)*2+jb(jk)

return

end subroutine gugadrt_check_rcas3
