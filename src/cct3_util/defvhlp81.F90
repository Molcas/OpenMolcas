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

subroutine defvhlp81(r2,v,dimr2b,dimr2a,dimr2c,dimva,dimvb,dimvc,adda,addc)
! this routine does
! V(a,b,c)aba = - R2(b,a,c)
! for syma>symc
!
! r2     - r2 matrix (I)
! v      - v matrix (O)
! dimr2b - dimension of b in R2 (I)
! dimr2a - dimension of a in R2 (I)
! dimr2c - dimension of c in R2 (I)
! dimva  - dimension of a in V (I)
! dimvb  - dimension of b in V (I)
! dimvc  - dimension of c in V (I)
! adda   - additional constant to a (I)
! addc   - additional constant to c (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimr2b, dimr2a, dimr2c, dimva, dimvb, dimvc, adda, addc
real(kind=wp), intent(in) :: r2(dimr2b,dimr2a,dimr2c)
real(kind=wp), intent(out) :: v(dimva,dimvb,dimvc)
integer(kind=iwp) :: a, ar2, c, cr2

do c=1,dimvc
  cr2 = c+addc
  do a=1,dimva
    ar2 = a+adda
    v(a,:,c) = -r2(1:dimvb,ar2,cr2)
  end do
end do

return

end subroutine defvhlp81
