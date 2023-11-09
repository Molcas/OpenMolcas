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

subroutine asonc10_cvb( &
#                      define _CALLING_
#                      include "ddasonc_interface.fh"
                      )

use casvb_global, only: cpu0, ipp10, iter10
use Definitions, only: wp, iwp, u6

implicit none
#include "ddasonc_interface.fh"
integer(kind=iwp) :: ivec
real(kind=wp), external :: tim_cvb

#include "macros.fh"
unused_var(sxc)

iter10 = iter10+1
if (ipp10 >= 2) then
  write(u6,'(/,a,i5,a,f10.3,a)') ' Davidson iteration',iter10,' at',tim_cvb(cpu0),' CPU seconds'
  write(u6,'(a)') ' -----------------------------------------------'
end if

axc(:,:) = c(:,:)
do ivec=1,nvec
  call hess_cvb(axc(:,ivec))
  call ddproj_cvb(axc(:,ivec),nprm)
end do

return

end subroutine asonc10_cvb
