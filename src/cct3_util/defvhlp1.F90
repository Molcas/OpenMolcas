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

subroutine defvhlp1(r1,v,dimr1a,dimr1bc,dimvab,dimvc,add)
! this routine does
! V(ab,c)xxx = R1(a,bc)-R1(b,ac) x=a,b
! for syma=symb=symc
!
! r1      - r1 matrix (I)
! v       - v matrix (O)
! dimr1a  - dimension of a (b,c) in R1 (I)
! dimr1bc - dimension of bc (ac) in R1 (I)
! dimvab  - dimension of ab in V (I)
! dimvc   - dimension of c (a,b) in V (I)
! add     - additional constant (I)
!           (# of singly occ in syma for alpha, 0 for beta)

use CCT3_global, only: nshf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimr1a, dimr1bc, dimvab, dimvc, add
real(kind=wp), intent(in) :: r1(dimr1a,dimr1bc)
real(kind=wp), intent(out) :: v(dimvab,dimvc)
integer(kind=iwp) :: a, ab, ab0, acr1, ar1, b, bk, c, cr1 !, bcr1

!do c=1,dimvc
!  cr1 = c+add
!  do a=2,dimvc
!    ar1 = a+add
!    ab0 = nshf(a)
!    do b=1,a-1
!      bcr1 = indab(b+add,cr1)
!      v(ab0+b,c) = r1(ar1,bcr1)
!    end do
!  end do
!end do

do c=1,dimvc
  cr1 = c+add
  do a=2,dimvc
    ar1 = a+add
    ab0 = nshf(a)
    if (c <= (a-2)) then
      ! b1 <= c1
      bk = cr1*(cr1-1)/2+add
      ! bcr1=cr1*(cr1-1)/2+add+b
      v(ab0+1:ab0+c,c) = r1(ar1,bk+1:bk+c)
      ! b1 > c1
      bk = (add+c+1)*(add+c)/2+cr1
      do b=c+1,a-1
        ! bcr1=(add+b)*(add+b-1)/2+cr1
        v(ab0+b,c) = r1(ar1,bk)
        bk = bk+add+b
      end do
    else
      bk = cr1*(cr1-1)/2+add
      ! b1 <= c1
      ! bcr1=cr1*(cr1-1)/2+add+b
      v(ab0+1:ab0+a-1,c) = r1(ar1,bk+1:bk+a-1)
    end if
  end do
end do

do c=1,dimvc
  cr1 = c+add
  do a=2,dimvc
    ab0 = nshf(a)
    ! acr1=indab(a+add,cr1)
    if ((a+add) > cr1) then
      acr1 = (a+add)*(a+add-1)/2+cr1
    else
      acr1 = cr1*(cr1-1)/2+a+add
    end if
    do b=1,a-1
      ab = ab0+b
      v(ab,c) = v(ab,c)-r1(b+add,acr1)
    end do
  end do
end do

return

end subroutine defvhlp1
