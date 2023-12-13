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
! Copyright (C) 2002, Roland Lindh                                     *
!***********************************************************************

subroutine Assemble_dVdB(NAInt,EFInt,nZeta,la,lb,A,B,C)
!***********************************************************************
!                                                                      *
!     Object: to assemble the derivative of the nuclear attraction     *
!             integrals with respect to the magnetic field.            *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemical Physics,                 *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb
real(kind=wp), intent(in) :: NAInt(nZeta*nTri_Elem1(la)*nTri_Elem1(lb)), A(3), B(3), C(3)
real(kind=wp), intent(inout) :: EFInt(nZeta*nTri_Elem1(la)*nTri_Elem1(lb),3)
integer(kind=iwp) :: iVec, nVec
real(kind=wp) :: EFInt_x, EFInt_y, EFInt_z, RAB(3)

RAB(:) = A-B

! Recombine in place!

nVec = size(EFInt,1)
do iVec=1,nVec
  EFInt_x = EFInt(iVec,1)
  EFInt_y = EFInt(iVec,2)
  EFInt_z = EFInt(iVec,3)
  EFInt(iVec,1) = RAB(2)*(EFInt_z+C(3)*NAInt(iVec))-RAB(3)*(EFInt_y+C(2)*NAInt(iVec))
  EFInt(iVec,2) = RAB(3)*(EFInt_x+C(1)*NAInt(iVec))-RAB(1)*(EFInt_z+C(3)*NAInt(iVec))
  EFInt(iVec,3) = RAB(1)*(EFInt_y+C(2)*NAInt(iVec))-RAB(2)*(EFInt_x+C(1)*NAInt(iVec))
end do

return

end subroutine Assemble_dVdB
