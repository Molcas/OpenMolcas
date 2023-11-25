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

subroutine asonc7_cvb( &
#                     define _CALLING_
#                     include "ddasonc_interface.fh"
                     )

use casvb_global, only: cpu0, ipp7, iter7, ograd
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
#include "ddasonc_interface.fh"
integer(kind=iwp) :: ivec
real(kind=wp), parameter :: thresh = 1.0e-15_wp
real(kind=wp), external :: ddot_, dnrm2_, tim_cvb

#include "macros.fh"
unused_var(sxc)

iter7 = iter7+1
if (ipp7 >= 2) then
  write(u6,'(/,a,i5,a,f10.3,a)') ' Davidson iteration',iter7,' at',tim_cvb(cpu0),' CPU seconds'
  write(u6,'(a)') ' -----------------------------------------------'
end if

do ivec=1,nvec
  axc(1,ivec) = ddot_(nprm-1,ograd,1,c(2:,ivec),1)
  axc(2:,ivec) = c(2:,ivec)
  ! Save Hessian application (& DNRM2 call) whenever possible:
  ! (C assumed to be normalized)
  if (abs(abs(c(1,ivec))-One) > thresh) then
    call hess_cvb(axc(2:,ivec))
  else if (dnrm2_(nprm-1,axc(2:,ivec),1) > thresh) then
    call hess_cvb(axc(2:,ivec))
  end if
  axc(2:,ivec) = axc(2:,ivec)+c(1,ivec)*ograd(1:nprm-1)
  call ddproj_cvb(axc(2:,ivec),nprm-1)
end do

return

end subroutine asonc7_cvb
