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

function ifns_cvb(nel1,nalf1,kbasis1)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ifns_cvb
integer(kind=iwp), intent(in) :: nel1, kbasis1
integer(kind=iwp), intent(inout) :: nalf1
integer(kind=iwp) :: ifn, iretval1, iretval2, nbet1, nsw

nbet1 = nel1-nalf1
if (nbet1 > nalf1) then
  nsw = nalf1
  nalf1 = nbet1
  nbet1 = nsw
end if
if (kbasis1 /= 6) then
  call icomb_cvb(nel1,nbet1,iretval1)
  call icomb_cvb(nel1,nbet1-1,iretval2)
  ifn = iretval1-iretval2
else
  call icomb_cvb(nel1,nalf1,ifn)
  if (nalf1 == nbet1) ifn = (ifn+1)/2
end if
ifns_cvb = ifn

return

end function ifns_cvb
