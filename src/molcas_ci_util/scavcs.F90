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

subroutine SCAVCS(VECO,VECI,INDEX,NDIM)
! Scatter vector with sign change
!
! vecO(abs(index(i))) = veci(i)*sign(index(i))

implicit real*8(A-H,O-Z)
dimension VECI(*), VECO(*), index(*)
intrinsic SIGN

do I=1,NDIM
  VECO(abs(index(I))) = VECI(I)*dble(sign(1,index(I)))
end do

return

end subroutine SCAVCS
