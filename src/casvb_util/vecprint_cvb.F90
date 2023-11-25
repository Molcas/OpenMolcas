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

subroutine vecprint_cvb(a,n)
! Prints vector A

use casvb_global, only: formMXP5, formMXP6, iprec, iwidth
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: a(n)
integer(kind=iwp) :: ibegin, iend, iform, nbuf
integer(kind=iwp), parameter :: mxbuf = 8

nbuf = min((iwidth-4)/(iprec+4),mxbuf)
if (nbuf == 7) nbuf = 6
iform = 0
do ibegin=1,n,nbuf
  iend = min(ibegin+nbuf-1,n)
  if (iform == 0) then
    write(u6,formMXP5) a(ibegin:iend)
  else
    write(u6,formMXP6) a(ibegin:iend)
  end if
end do

return

end subroutine vecprint_cvb
