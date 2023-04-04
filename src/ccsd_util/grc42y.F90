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

subroutine grc42y(mapda,mapdb,mapdc,mapia,mapib,mapic,mvec,ssa,ssb,possc0,ix)

use ccsd_global, only: dimm, mmul, nsym
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: mapda(0:512,6), mapdb(0:512,6), mapdc(0:512,6), mapia(8,8,8), mapib(8,8,8), mapic(8,8,8), mvec(4096,7), ssa, &
                     ssb, possc0, ix
integer(kind=iwp) :: ia, ib, iy, nhelp1, nhelp2, nhelp21, nhelp22, nhelp4, nhelp41, nhelp42, ntest1, ntest2, possct, sa1, sa134, &
                     sa2, sa3, sa34, sa4, sb1, sb2

! structure A(pq,rs)*B(rs)=YC(pq)

!1.0 prepare mapdc,mapic

if ((mapda(0,6) == 1) .or. (mapda(0,6) == 4)) then
  ntest1 = 1
else
  ntest1 = 0
end if

call grc0(2,ntest1,mapda(0,1),mapda(0,2),0,0,mmul(ssa,ssb),possc0,possct,mapdc,mapic)

!1.1 define limitations - p>q,r,s must be tested - ntest1
!                       - p,q,r>s must be tested - ntest2

if ((mapda(0,6) == 1) .or. (mapda(0,6) == 4)) then
  ntest1 = 1
else
  ntest1 = 0
end if

if ((mapda(0,6) == 3) .or. (mapda(0,6) == 4)) then
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

    !1.3 def mvec,mapdc and mapdi

    ia = mapia(sa1,sa2,sa3)
    ib = mapib(sb1,1,1)
    iy = mapic(sa1,1,1)

    ! yes/no
    if ((mapda(ia,2) <= 0) .or. (mapdb(ib,2) <= 0)) cycle
    nhelp1 = 1

    ! rowA
    nhelp21 = dimm(mapda(0,1),sa1)
    nhelp22 = dimm(mapda(0,2),sa2)
    if ((ntest1 == 1) .and. (sa1 == sa2)) then
      nhelp2 = nhelp21*(nhelp21-1)/2
    else
      nhelp2 = nhelp21*nhelp22
    end if

    ! sum
    nhelp41 = dimm(mapda(0,3),sa3)
    nhelp42 = dimm(mapda(0,4),sa4)
    if ((ntest2 == 1) .and. (sa3 == sa4)) then
      nhelp4 = nhelp41*(nhelp41-1)/2
    else
      nhelp4 = nhelp41*nhelp42
    end if

    mvec(ix,1) = nhelp1
    mvec(ix,2) = mapda(ia,1)
    mvec(ix,3) = mapdb(ib,1)
    mvec(ix,4) = mapdc(iy,1)
    mvec(ix,5) = nhelp2
    mvec(ix,6) = nhelp4
    mvec(ix,7) = 0

    ix = ix+1

  end do
end do
ix = ix-1

return

end subroutine grc42y
