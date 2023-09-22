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

subroutine asonc10_cvb(c,axc,dum1,nvec,nprm)

use casvb_global, only: ipp, iter
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nvec, nprm
real(kind=wp) :: c(nprm,nvec), axc(nprm,nvec), dum1
#include "main_cvb.fh"
integer(kind=iwp) :: ivec
real(kind=wp), external :: tim_cvb

iter = iter+1
if (ipp >= 2) then
  write(u6,'(/,a,i5,a,f10.3,a)') ' Davidson iteration',iter,' at',tim_cvb(cpu0),' CPU seconds'
  write(u6,'(a)') ' -----------------------------------------------'
end if

do ivec=1,nvec
  call fmove_cvb(c(1,ivec),axc(1,ivec),nprm)
  call hess_cvb(axc(1,ivec))
  call ddproj_cvb(axc(1,ivec),nprm)
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_real(dum1)

end subroutine asonc10_cvb
