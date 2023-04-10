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

implicit none
integer ncLoc, nc, ncOff, dim_1
real*8 Ll(1:ncLoc,1:dim_1)
real*8 Lg(1:nc,1:dim_1)
! help variables
integer m, p

do p=1,dim_1
  do m=1,ncLoc
    Lg(ncOff+m,p) = Ll(m,p)
  end do
end do

return

end subroutine InsL
