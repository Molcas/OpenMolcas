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

subroutine bikset_cvb(aikcof,bikcof,nel,nalf,i2s,ndet,ifns,kbasis,share,iprint)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nel, nalf, i2s, ndet, ifns, kbasis, iprint
real(kind=wp), intent(inout) :: aikcof(ndet,ifns)
real(kind=wp), intent(out) :: bikcof(ndet,ifns)
logical(kind=iwp), intent(in) :: share
integer(kind=iwp) :: nalf_use, ndet_use
real(kind=wp), allocatable :: atmp(:,:), btmp(:,:)

if (i2s /= 2*nalf-nel) then
  nalf_use = (i2s+nel)/2
  call icomb_cvb(nel,nalf_use,ndet_use)
  call mma_allocate(atmp,ndet_use,ifns,label='atmp')
  call mma_allocate(btmp,ndet_use,ifns,label='btmp')
  call biksmain_cvb(atmp,btmp,nel,nalf_use,ndet_use,ifns,kbasis,.false.,iprint)
  call sminus_cvb(atmp,bikcof,nel,nalf_use,nalf,ifns)
  if (.not. share) call sminus_cvb(btmp,aikcof,nel,nalf_use,nalf,ifns)
  call mma_deallocate(atmp)
  call mma_deallocate(btmp)
else
  call biksmain_cvb(aikcof,bikcof,nel,nalf,ndet,ifns,kbasis,share,iprint)
end if

return

end subroutine bikset_cvb
