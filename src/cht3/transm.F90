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

subroutine TRANSM(V1,V2,ID1,ID2)

implicit none
integer ID1, ID2, I
real*8 V1, V2, ONE
parameter(ONE=1.d0)
dimension V1(ID1,ID2), V2(ID2,ID1)

! usual transposition
if ((ID1 == 0) .or. (ID2 == 0)) return
do I=1,ID1
  call DCOPY_(ID2,V1(I,1),ID1,V2(1,I),1)
end do

return

entry TRANSM_A(V1,V2,ID1,ID2)
! transposition of the matrix + add to the target matrix

if ((ID1 == 0) .or. (ID2 == 0)) return
do I=1,ID1
  call DAXPY_(ID2,ONE,V1(I,1),ID1,V2(1,I),1)
end do

return

end subroutine TRANSM
