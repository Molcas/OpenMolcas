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

! Rotate multipole.
subroutine Rotation_qmstat(iL,dMul,Rotte,Sigge)

use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

! Maximum multipole implemented
#define _MxM_ 2

implicit none
integer(kind=iwp), intent(in) :: iL
real(kind=wp), intent(inout) :: dMul(nTri_Elem1(_MxM_))
real(kind=wp), intent(in) :: Rotte(3,3), Sigge
real(kind=wp) :: d1, d2, d3, dMTrans(6), Sig, TD(6,6)
integer(kind=iwp) :: i, j
#include "warnings.h"

if (iL == 0) then
  ! Charge, trivial to rotate.

else if (iL == 1) then
  ! Dipole, transforms as a vector. Sigge controls that if the
  ! multipole is located not in origin, but at the other end,
  ! i.e. molecule A, then any odd occurrence of z should be
  ! mirrored. Applies for the quadrupole as well, see below.

  d1 = dMul(1)
  d2 = dMul(2)
  d3 = dMul(3)
  dMul(:) = Rotte(:,1)*d1+Rotte(:,2)*d2+Rotte(:,3)*d3
  dMul(3) = Sigge*dMul(3)
else if (iL == 2) then
  ! Quadrupole, transforms as a quadratic form. Also, transform to spherical representation.

  ! Compute the transformation matrix for second-moments.

  call M2Trans(Rotte,TD)

  ! Transform. Sigge is explained above.

  dMTrans(:) = Zero
  do i=1,6
    do j=1,6
      dMTrans(i) = dMTrans(i)+TD(i,j)*dMul(j)
    end do
  end do
  do i=1,6
    Sig = One
    if ((i == 3) .or. (i == 5)) Sig = Sigge
    dMul(i) = dMTrans(i)*Sig
  end do

  ! Go to spherical representation.

  call Spherical(dMul)
else
  write(u6,*) 'Nope!, Error in sl_grad'
  call Quit(_RC_IO_ERROR_READ_)
end if

return

end subroutine Rotation_qmstat
