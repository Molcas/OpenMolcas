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

real*8 function INPROD(A,B,NDIM)
! CALCULATE SCALAR PRODUCT BETWEEN TWO VECTORS A,B

implicit real*8(A-H,O-Z)
dimension A(*), B(*)

INPROD = 0.0d0
do I=1,NDIM
  INPROD = INPROD+A(I)*B(I)
end do

end function INPROD
