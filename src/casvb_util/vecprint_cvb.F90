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

implicit real*8(a-h,o-z)
#include "print_cvb.fh"
#include "formats_cvb.fh"
parameter(mxbuf=8)
dimension a(n)

nbuf = min((iwidth-4)/(iprec+4),mxbuf)
if (nbuf == 7) nbuf = 6
iform = 0
do ibegin=1,n,nbuf
  iend = min(ibegin+nbuf-1,n)
  if (iform == 0) then
    write(6,formMXP5) (a(i),i=ibegin,iend)
  else
    write(6,formMXP6) (a(i),i=ibegin,iend)
  end if
end do

return

end subroutine vecprint_cvb
