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

subroutine defvhlp3(r1,r2,v,dimr1a,dimr1b,dimr1c,dimr2ac,dimva,dimvb,dimvc,adda,addb,addc)
! this routine does
! V(a,b,c)xxx = R1(a,c,b)-R2(b,ac) x=a,b
! for syma/=symb symc==syma
!
! r1      - r1 matrix (I)
! r2      - r2 matrix (I)
! v       - v matrix (O)
! dimr1a  - dimension of a in R1 (I)
! dimr1b  - dimension of b in R1 (I)
! dimr1c  - dimension of c in R1 (I)
! dimr2ac - dimension of ac in R2 (I)
! dimva   - dimension of a in V (I)
! dimvb   - dimension of b in V (I)
! dimvc   - dimension of c in V (I)
! adda    - additional constant to a (I)
! addb    - additional constant to b (I)
! addc    - additional constant to c (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimr1a, dimr1b, dimr1c, dimr2ac, dimva, dimvb, dimvc, adda, addb, addc
real(kind=wp), intent(in) :: r1(dimr1a,dimr1c,dimr1b), r2(dimr1b,dimr2ac)
real(kind=wp), intent(out) :: v(dimva,dimvb,dimvc)
integer(kind=iwp) :: a, acr2, b, br1, c, cr1, cr2

do b=1,dimvb
  br1 = b+addb
  do c=1,dimvc
    cr1 = c+addc
    v(:,b,c) = r1(adda+1:adda+dimva,cr1,br1)
  end do
end do

do c=1,dimvc
  cr2 = c+addc
  do a=1,dimvc
    !acr2 = indab(a+adda,cr2)
    if ((a+adda) >= cr2) then
      acr2 = (a+adda)*(a+adda-1)/2+cr2
    else
      acr2 = cr2*(cr2-1)/2+a+adda
    end if
    v(a,:,c) = v(a,:,c)-r2(addb+1:addb+dimvb,acr2)
  end do
end do

return

end subroutine defvhlp3
