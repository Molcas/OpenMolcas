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

subroutine party2_cvb(iperm,n,party)

implicit real*8(a-h,o-z)
dimension iperm(n)

! Following caters for non-contiguous integers:
ntransp = 0
100 continue
do i=1,n-1
  if (iperm(i) > iperm(i+1)) then
    ntransp = ntransp+1
    iswp = iperm(i)
    iperm(i) = iperm(i+1)
    iperm(i+1) = iswp
    do j=i-1,1,-1
      if (iperm(j) > iperm(j+1)) then
        ntransp = ntransp+1
        iswp = iperm(j)
        iperm(j) = iperm(j+1)
        iperm(j+1) = iswp
      end if
    end do
    goto 100
  end if
end do
if (mod(ntransp-n,2) == 0) then
  party = 1d0
else
  party = -1d0
end if

return

end subroutine party2_cvb
