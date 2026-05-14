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

subroutine free2all_cvb(vecfrom,vecto,nvec)

use casvb_global, only: nfr, nfrorb, npr, nprorb, nprvb, orbfr_is_unit, trprm
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nvec
real(kind=wp), intent(in) :: vecfrom(nfr,nvec)
real(kind=wp), intent(out) :: vecto(npr,nvec)
integer(kind=iwp) :: ivec

do ivec=1,nvec
  if (.not. orbfr_is_unit) then
    call mxatb_cvb(trprm,vecfrom(:,ivec),nprorb,nfrorb,1,vecto(:,ivec))
  else
    if (nprorb > 0) vecto(1:nprorb,ivec) = vecfrom(1:nprorb,ivec)
  end if
  if (nprvb > 0) vecto(nprorb+1:nprorb+nprvb,ivec) = vecfrom(nfrorb+1:nfrorb+nprvb,ivec)
end do

return

end subroutine free2all_cvb
