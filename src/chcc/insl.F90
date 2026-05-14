!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine InsL(Ll,Lg,ncLoc,nc,ncOff,dim_1)
! this routine does:
! insert Llocal(ml,dim_1) into Lglobal(m,dim_1)
! on a corresponding place

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ncLoc, nc, ncOff, dim_1
real(kind=wp), intent(in) :: Ll(ncLoc,dim_1)
real(kind=wp), intent(inout) :: Lg(nc,dim_1)

Lg(ncOff+1:ncOff+ncLoc,:) = Ll(:,:)

return

end subroutine InsL
