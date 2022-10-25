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

subroutine defvhlp61(r1,v,dimr1a,dimr1b,dimr1c,dimva,dimvb,dimvc,adda)
! this routine does
! V(a,b,c)abb = R1(a,b,c)
! for symb>symc
!
! r1     - r1 matrix (I)
! v      - v matrix (O)
! dimr1a - dimension of a in R1 (I)
! dimr1b - dimension of b in R1 (I)
! dimr1c - dimension of c in R1 (I)
! dimva  - dimension of a in V (I)
! dimvb  - dimension of b in V (I)
! dimvc  - dimension of c in V (I)
! adda   - additional constant to a (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimr1a, dimr1b, dimr1c, dimva, dimvb, dimvc, adda
real(kind=wp), intent(in) :: r1(dimr1a,dimr1b,dimr1c)
real(kind=wp), intent(out) :: v(dimva,dimvb,dimvc)

v(:,:,:) = r1(adda+1:adda+dimva,1:dimvb,1:dimvc)

return

end subroutine defvhlp61
