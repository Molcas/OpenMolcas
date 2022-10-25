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

subroutine defvhlp52(r1,r2,v,dimr1a,dimr1b,dimr1c,dimva,dimvb,dimvc,adda,addb,addc)
! this routine does
! V(a,b,c)xxx = R1(a,b,c)-R2(b,a,c) x=a,b
! for syma>symc, symb<symc, (syma>symb)
!
! r1     - r1 matrix (I)
! r2     - r2 matrix (I)
! v      - v matrix (O)
! dimr1a - dimension of a in R (I)
! dimr1b - dimension of b in R (I)
! dimr1c - dimension of c in R (I)
! dimva  - dimension of a in V (I)
! dimvb  - dimension of b in V (I)
! dimvc  - dimension of c in V (I)
! adda   - additional constant to a (I)
! addb   - additional constant to b (I)
! addc   - additional constant to c (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimr1a, dimr1b, dimr1c, dimva, dimvb, dimvc, adda, addb, addc
real(kind=wp), intent(in) :: r1(dimr1a,dimr1c,dimr1b), r2(dimr1b,dimr1a,dimr1c)
real(kind=wp), intent(out) :: v(dimva,dimvb,dimvc)
integer(kind=iwp) :: b, br1, br2, c, cr1, cr2

do b=1,dimvb
  br1 = b+addb
  do c=1,dimvc
    cr1 = c+addc
    v(:,b,c) = r1(adda+1:adda+dimva,cr1,br1)
  end do
end do

do c=1,dimvc
  cr2 = c+addc
  do b=1,dimvb
    br2 = b+addb
    v(:,b,c) = v(:,b,c)-r2(br2,adda+1:adda+dimva,cr2)
  end do
end do

return

end subroutine defvhlp52
