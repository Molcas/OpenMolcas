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

subroutine CHO_INVPCK(IJ,I,J,LOW)
!
! Purpose: invert triangular packing index ij to
!          rectangular indices i and j.
!          Flag LOW specifies packing convention:
!          LOW = T: i>=j
!          LOW = F: i<=j

implicit real*8(a-h,o-z)
logical LOW
parameter(ONE=1.0d0,TWO=2.0d0,THREE=3.0d0,EIGHT=8.0d0)

if (IJ > 0) then

  XX = EIGHT*dble(IJ)-THREE
  XI = (ONE+sqrt(XX))/TWO
  I = int(XI)
  J = IJ-I*(I-1)/2

  if (.not. LOW) then
    ITMP = I
    I = J
    J = ITMP
  end if

else

  I = -1
  J = -2

end if

end subroutine CHO_INVPCK
