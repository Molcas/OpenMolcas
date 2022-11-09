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

subroutine JACORD(HH,EIGVEC,NVEC,NDIM)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NVEC, NDIM
real(kind=wp), intent(inout) :: HH(*), EIGVEC(NDIM,NVEC)
integer(kind=iwp) :: I, II, IMIN, J, JJ, K
real(kind=wp) :: EI, EJ, EMIN, SWAP
real(kind=wp), parameter :: ThrZ = 1.0e-14_wp

do I=1,NVEC-1
  II = nTri_Elem(I)
  EI = HH(II)
  EMIN = EI
  IMIN = I
  do J=I+1,NVEC
    JJ = nTri_Elem(J)
    EJ = HH(JJ)
    if ((EJ >= EMIN) .or. (abs(EJ-EMIN) < ThrZ)) cycle
    EMIN = EJ
    IMIN = J
  end do
  if (IMIN == I) cycle
  HH(II) = EMIN
  JJ = nTri_Elem(IMIN)
  HH(JJ) = EI
  do K=1,NDIM
    SWAP = EIGVEC(K,I)
    EIGVEC(K,I) = EIGVEC(K,IMIN)
    EIGVEC(K,IMIN) = SWAP
  end do
end do

return

end subroutine JACORD
