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

implicit real*8(A-H,O-Z)
dimension HH(*), EIGVEC(NDIM,NVEC)

ThrZ = 1.0D-14
do I=1,NVEC-1
  II = (I*(I+1))/2
  EI = HH(II)
  EMIN = EI
  IMIN = I
  do J=I+1,NVEC
    JJ = (J*(J+1))/2
    EJ = HH(JJ)
    if ((EJ >= EMIN) .or. (abs(EJ-EMIN) < ThrZ)) cycle
    EMIN = EJ
    IMIN = J
  end do
  if (IMIN == I) cycle
  HH(II) = EMIN
  JJ = (IMIN*(IMIN+1))/2
  HH(JJ) = EI
  do K=1,NDIM
    SWAP = EIGVEC(K,I)
    EIGVEC(K,I) = EIGVEC(K,IMIN)
    EIGVEC(K,IMIN) = SWAP
  end do
end do

return

end subroutine JACORD
