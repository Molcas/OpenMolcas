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

function ioemrg2_cvb(ia1,na1,ia2,na2,ia12)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ioemrg2_cvb
integer(kind=iwp) :: na1, ia1(na1), na2, ia2(na2), ia12(*)
integer(kind=iwp) :: i, ioe, n1, n12, n2

n1 = 1
n2 = 1
ioe = 0
n12 = 1
do
  if (n1 > na1) then
    ioemrg2_cvb = 1-2*mod(ioe,2)
    do i=n2,na2
      ia12(n12) = ia2(i)
      n12 = n12+1
    end do
    exit
  else if (n2 > na2) then
    ioe = ioe+(na1-n1+1)*na2
    ioemrg2_cvb = 1-2*mod(ioe,2)
    do i=n1,na1
      ia12(n12) = ia1(i)
      n12 = n12+1
    end do
    exit
  end if
  if (ia1(n1) < ia2(n2)) then
    ioe = ioe+n2-1
    ia12(n12) = ia1(n1)
    n1 = n1+1
    n12 = n12+1
  else if (ia1(n1) > ia2(n2)) then
    ia12(n12) = ia2(n1)
    n2 = n2+1
    n12 = n12+1
  else if (ia1(n1) == ia2(n2)) then
    ioemrg2_cvb = 0
    exit
  end if
end do

return

end function ioemrg2_cvb
