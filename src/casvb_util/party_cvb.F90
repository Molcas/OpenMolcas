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

function party_cvb(iperm,n)
! Returns parity of permutation

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: party_cvb
integer(kind=iwp), intent(in) :: n, iperm(n)
integer(kind=iwp) :: i, iswp, j, ntransp
integer(kind=iwp), allocatable :: tmp(:)
logical(kind=iwp) :: done

call mma_allocate(tmp,n,label='tmp')
tmp(:) = iperm(:)

! Following caters for non-contiguous integers:
ntransp = 0
do
  done = .false.
  do i=1,n-1
    if (tmp(i) > tmp(i+1)) then
      ntransp = ntransp+1
      iswp = tmp(i)
      tmp(i) = tmp(i+1)
      tmp(i+1) = iswp
      do j=i-1,1,-1
        if (tmp(j) > tmp(j+1)) then
          ntransp = ntransp+1
          iswp = tmp(j)
          tmp(j) = tmp(j+1)
          tmp(j+1) = iswp
        end if
      end do
      done = .true.
      exit
    end if
  end do
  if (.not. done) exit
end do
if (mod(ntransp-n,2) == 0) then
  party_cvb = One
else
  party_cvb = -One
end if

call mma_deallocate(tmp)

return

end function party_cvb
