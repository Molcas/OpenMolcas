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

subroutine GATVCS(VECO,VECI,INDEX,NDIM)
! Gather vector allowing for sign change
!
! VECO(I) = VECI(INDEX(I))

implicit real*8(A-H,O-Z)
dimension VECI(*), VECO(*), index(*)
intrinsic SIGN

do I=1,NDIM
  VECO(I) = VECI(abs(index(I)))*dble(sign(1,index(I)))
end do

return

end subroutine GATVCS
