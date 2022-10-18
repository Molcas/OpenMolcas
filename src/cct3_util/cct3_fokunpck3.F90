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

subroutine cct3_fokunpck3(fok,fai,dimfok,dimfa,dimfi)
! this routine distributes (Fok - dp) to Fai
! fok    - Fok matrix (I)
! fai    - Fai matrix (O)
! dimfok - dimension for Fok matrix - norb (I)
! dimfa  - dimension of virtuals - nv (I)
! dimfi  - dimension of occupied - no (I)

integer dimfok, dimfa, dimfi

real*8 fok(1:dimfok,1:dimfok)
real*8 fai(1:dimfa,1:dimfi)
! help variables
integer a, i

!1 distribute Fok to Fai
do i=1,dimfi
  do a=1,dimfa
    fai(a,i) = fok(dimfi+a,i)
  end do
end do

return

end subroutine cct3_fokunpck3
