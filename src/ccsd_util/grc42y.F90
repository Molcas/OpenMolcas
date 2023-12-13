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

subroutine grc42y(a,b,c,mvec,ssa,ssb,ix)

use ccsd_global, only: dimm, Map_Type, mmul, nsym
use Definitions, only: iwp

implicit none
type(Map_Type), intent(in) :: a, b
type(Map_Type), intent(inout) :: c
integer(kind=iwp), intent(out) :: mvec(4096,7), ix
integer(kind=iwp), intent(in) :: ssa, ssb
integer(kind=iwp) :: ia, ib, iy, nhelp1, nhelp2, nhelp21, nhelp22, nhelp4, nhelp41, nhelp42, ntest1, ntest2, posct, sa1, sa134, &
                     sa2, sa3, sa34, sa4, sb1, sb2

! structure A(pq,rs)*B(rs)=YC(pq)

!1.0 prepare c%d,c%i

if ((a%d(0,6) == 1) .or. (a%d(0,6) == 4)) then
  ntest1 = 1
else
  ntest1 = 0
end if

call grc0(2,ntest1,a%d(0,1),a%d(0,2),0,0,mmul(ssa,ssb),posct,c)

!1.1 define limitations - p>q,r,s must be tested - ntest1
!                       - p,q,r>s must be tested - ntest2

if ((a%d(0,6) == 1) .or. (a%d(0,6) == 4)) then
  ntest1 = 1
else
  ntest1 = 0
end if

if ((a%d(0,6) == 3) .or. (a%d(0,6) == 4)) then
  ntest2 = 1
else
  ntest2 = 0
end if

!1.2 def symm states and test the limitations

ix = 1
do sb1=1,nsym
  sa3 = sb1

  sb2 = mmul(ssb,sb1)
  sa4 = sb2
  sa34 = mmul(sa3,sa4)
  ! Meggie out
  if ((ntest2 == 1) .and. (sb1 < sb2)) cycle

  do sa1=1,nsym
    sa134 = mmul(sa1,sa34)

    sa2 = mmul(ssa,sa134)
    ! Meggie out
    if ((ntest1 == 1) .and. (sa1 < sa2)) cycle

    !1.3 def mvec,c%d and c%i

    ia = a%i(sa1,sa2,sa3)
    ib = b%i(sb1,1,1)
    iy = c%i(sa1,1,1)

    ! yes/no
    if ((a%d(ia,2) <= 0) .or. (b%d(ib,2) <= 0)) cycle
    nhelp1 = 1

    ! rowA
    nhelp21 = dimm(a%d(0,1),sa1)
    nhelp22 = dimm(a%d(0,2),sa2)
    if ((ntest1 == 1) .and. (sa1 == sa2)) then
      nhelp2 = nhelp21*(nhelp21-1)/2
    else
      nhelp2 = nhelp21*nhelp22
    end if

    ! sum
    nhelp41 = dimm(a%d(0,3),sa3)
    nhelp42 = dimm(a%d(0,4),sa4)
    if ((ntest2 == 1) .and. (sa3 == sa4)) then
      nhelp4 = nhelp41*(nhelp41-1)/2
    else
      nhelp4 = nhelp41*nhelp42
    end if

    mvec(ix,1) = nhelp1
    mvec(ix,2) = a%d(ia,1)
    mvec(ix,3) = b%d(ib,1)
    mvec(ix,4) = c%d(iy,1)
    mvec(ix,5) = nhelp2
    mvec(ix,6) = nhelp4
    mvec(ix,7) = 0

    ix = ix+1

  end do
end do
ix = ix-1

return

end subroutine grc42y
