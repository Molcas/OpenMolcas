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

subroutine defvhlp22(r1,v,dimr1a,dimr1c,dimvab,dimva,dimvc,adda,addc)
! this routine does
! V(ab,c)xxx = R1(a,c,b)-R1(b,c,a) x=a,b
! for syma=symb<symc
!
! r1     - r1 matrix (I)
! v      - v matrix (O)
! dimr1a - dimension of a (b) in R1 (I)
! dimr1c - dimension of c in R1 (I)
! dimvab - dimension of ab in V (I)
! dimva  - dimension of a (b) in V (I)
! dimvc  - dimension of c in V (I)
! adda   - additional constant to a (b) (I)
! addc   - additional constant to c (I)

use CCT3_global, only: nshf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimr1a, dimr1c, dimvab, dimva, dimvc, adda, addc
real(kind=wp), intent(in) :: r1(dimr1a,dimr1c,dimr1a)
real(kind=wp), intent(out) :: v(dimvab,dimvc)
integer(kind=iwp) :: a, ab, ab0, ar1, b, c, cr1

do c=1,dimvc
  cr1 = c+addc
  do a=2,dimva
    ar1 = a+adda
    ab0 = nshf(a)
    v(ab0+1:ab0+a-1,c) = r1(ar1,cr1,adda+1:adda+a-1)
  end do
end do

do a=2,dimva
  ar1 = a+adda
  ab0 = nshf(a)
  do c=1,dimvc
    cr1 = c+addc
    do b=1,a-1
      ab = ab0+b
      v(ab,c) = v(ab,c)-r1(b+adda,cr1,ar1)
    end do
  end do
end do

return

end subroutine defvhlp22
