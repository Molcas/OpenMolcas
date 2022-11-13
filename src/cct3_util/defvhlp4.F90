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

subroutine defvhlp4(r1,r2,v,dimr1a,dimr1bc,dimr2b,dimr2c,dimva,dimvb,dimvc,adda,addb,addc)
! this routine does
! V(a,b,c)xxx = R1(a,bc)-R2(b,a,c) x=a,b
! for syma/=symb symc==symb
!
! r1      - r1 matrix (I)
! r2      - r2 matrix (I)
! v       - v matrix (O)
! dimr1a  - dimension of a in R1 (I)
! dimr1bc - dimension of bc in R1 (I)
! dimr2b  - dimension of b in R2 (I)
! dimr2c  - dimension of c in R2 (I)
! dimva   - dimension of a in V (I)
! dimvb   - dimension of b in V (I)
! dimvc   - dimension of c in V (I)
! adda    - additional constant to a (I)
! addb    - additional constant to b (I)
! addc    - additional constant to c (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimr1a, dimr1bc, dimr2b, dimr2c, dimva, dimvb, dimvc, adda, addb, addc
real(kind=wp), intent(in) :: r1(dimr1a,dimr1bc), r2(dimr2b,dimr1a,dimr2c)
real(kind=wp), intent(out) :: v(dimva,dimvb,dimvc)
integer(kind=iwp) :: b, bcr1, br2, c, cr1, cr2

do c=1,dimvc
  cr1 = c+addc
  do b=1,dimvb
    !bcr1 = indab(b+addb,cr1)
    if ((b+addb) > cr1) then
      bcr1 = (b+addb)*(b+addb-1)/2+cr1
    else
      bcr1 = cr1*(cr1-1)/2+b+addb
    end if
    v(:,b,c) = r1(adda+1:adda+dimva,bcr1)
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

end subroutine defvhlp4
