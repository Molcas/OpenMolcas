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

subroutine aikcof_cvb(aikcof,bikcof,ndet,ifns,kbasis,share)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ndet, ifns, kbasis
real(kind=wp), intent(inout) :: aikcof(ndet,ifns)
real(kind=wp), intent(in) :: bikcof(ndet,ifns)
logical(kind=iwp), intent(in) :: share
real(kind=wp), allocatable :: sovr(:,:)

if (kbasis == 6) return

! Generate mapping from determinants to spin functions
! (If KBASIS<=2 then AIKCOF=BIKCOF and they (probably) share memory)
if (kbasis > 2) then
  call mma_allocate(sovr,ifns,ifns,label='sovr')
  call mxattb_cvb(bikcof,bikcof,ifns,ndet,ifns,sovr)
  call mxinv_cvb(sovr,ifns)
  call mxatb_cvb(bikcof,sovr,ndet,ifns,ifns,aikcof)
  call mma_deallocate(sovr)
else if (.not. share) then
  aikcof(:,:) = bikcof(:,:)
end if

return

end subroutine aikcof_cvb
