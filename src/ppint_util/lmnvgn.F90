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

subroutine lmnvgn(lmn1u,lmnv)
! Generate Cartesian Gaussian exponent array.
! lmnv(*,*) = exponents of the cartesian gaussian basis functions.
!             s   p   d   f   g   h   i
!       lmn = 0   1   2   3   4   5   6
!    numxyz = 1,  3,  6, 10, 15  21  28      = ((lmn+1)*(lmn+2))/2

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lmn1u
integer(kind=iwp), intent(out) :: lmnv(3,lmn1u*(lmn1u+1)*(lmn1u+2)/6)
integer(kind=iwp) :: ix, iy, iz, lmn, ndx

ndx = 0
do lmn=0,lmn1u-1
  do ix=lmn,0,-1
    do iy=lmn-ix,0,-1
      iz = lmn-ix-iy
      ndx = ndx+1
      lmnv(1,ndx) = ix
      lmnv(2,ndx) = iy
      lmnv(3,ndx) = iz
    end do
  end do
end do

return

end subroutine lmnvgn
