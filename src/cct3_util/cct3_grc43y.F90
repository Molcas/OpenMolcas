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

subroutine cct3_grc43y(mapda,mapdb,mapdc,mapia,mapib,mapic,mvec,ssa,ssb,posc0,ix)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: mapda(0:512,6), mapdb(0:512,6), mapdc(0:512,6), mapia(8,8,8), mapib(8,8,8), mapic(8,8,8), mvec(4096,7), ssa, &
                     ssb, posc0, ix
#include "t31.fh"
integer(kind=iwp) :: ia, ib, iy, nhelp1, nhelp2, nhelp4, nhelp41, nhelp42, nhelp43, nsymb2, ntest1, ntest2, posct, sa1, sa2, sa23, &
                     sa234, sa3, sa4, sb1, sb12, sb2, sb3

! structure A(p,qrs)*B(qrs)=YC(p)

!1.0 prepare mapdc,mapic

call cct3_grc0(1,0,mapda(0,1),0,0,0,mmul(ssa,ssb),posc0,posct,mapdc,mapic)

!1.1 define limitations - p,q>r,s must be tested - ntest1
!                       - p,q,r>s must be tested - ntest2

if (mapda(0,6) == 2) then
  ntest1 = 1
else
  ntest1 = 0
end if

if (mapda(0,6) == 3) then
  ntest2 = 1
else
  ntest2 = 0
end if

!1.2 def symm states and test the limitations

ix = 1
do sb1=1,nsym
  sa2 = sb1
  if (ntest1 == 1) then
    nsymb2 = sb1
  else
    nsymb2 = nsym
  end if

  do sb2=1,nsymb2
    sb12 = mmul(sb1,sb2)
    sa3 = sb2
    sa23 = mmul(sa2,sa3)

    sb3 = mmul(ssb,sb12)
    sa4 = sb3
    sa234 = mmul(sa23,sa4)
    ! Meggie out
    if ((ntest2 == 1) .and. (sb2 < sb3)) cycle

    sa1 = mmul(ssa,sa234)
    ! Meggie out
    if ((ntest1 == 1) .and. (sb1 < sb2)) cycle

    !1.3 def mvec,mapdc and mapdi

    ia = mapia(sa1,sa2,sa3)
    ib = mapib(sb1,sb2,1)
    iy = mapic(1,1,1)

    ! yes/no
    if ((mapda(ia,2) > 0) .and. (mapdb(ib,2) > 0)) then
      nhelp1 = 1
    else
      cycle
    end if

    ! rowA
    nhelp2 = dimm(mapda(0,1),sa1)

    ! sum
    nhelp41 = dimm(mapda(0,2),sa2)
    nhelp42 = dimm(mapda(0,3),sa3)
    nhelp43 = dimm(mapda(0,4),sa4)
    if ((ntest1 == 1) .and. (sa2 == sa3)) then
      nhelp4 = nhelp41*(nhelp41-1)*nhelp43/2
    else if ((ntest2 == 1) .and. (sa3 == sa4)) then
      nhelp4 = nhelp41*nhelp42*(nhelp42-1)/2
    else
      nhelp4 = nhelp41*nhelp42*nhelp43
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

end subroutine cct3_grc43y
