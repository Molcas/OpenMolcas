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

subroutine sortindxi_cvb(n,arrin,indx)
! -- purpose: indexing arrin so that arrin(indx(j)), j=1..n is in
!    ascending numerical order.
!    method is heapsort, see also subroutine hpsort.
! -- taken from numerical recipes, p 233.

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: n, arrin(n)
integer(kind=iwp), intent(out) :: indx(n)
integer(kind=iwp) :: i, indxt, ir, j, l, q

do j=1,n
  indx(j) = j
end do
if (n <= 1) return
l = n/2+1
ir = n
do
  if (l > 1) then
    l = l-1
    indxt = indx(l)
    q = arrin(indxt)
  else
    indxt = indx(ir)
    q = arrin(indxt)
    indx(ir) = indx(1)
    ir = ir-1
    if (ir == 1) then
      indx(1) = indxt
      exit
    end if
  end if
  i = l
  j = 2*l
  do while (j <= ir)
    if (j < ir) then
      if (arrin(indx(j)) < arrin(indx(j+1))) j = j+1
    end if
    if (q < arrin(indx(j))) then
      indx(i) = indx(j)
      i = j
      j = 2*j
    else
      j = ir+1
    end if
  end do
  indx(i) = indxt
end do

return

end subroutine sortindxi_cvb
