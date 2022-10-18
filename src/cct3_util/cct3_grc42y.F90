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

subroutine cct3_grc42y(mapda,mapdb,mapdc,mapia,mapib,mapic,mvec,ssa,ssb,possc0,ix)

#include "t31.fh"
integer mapda(0:512,1:6)
integer mapdb(0:512,1:6)
integer mapdc(0:512,1:6)
integer mapia(1:8,1:8,1:8)
integer mapib(1:8,1:8,1:8)
integer mapic(1:8,1:8,1:8)
integer mvec(1:4096,1:7)
integer possc0
integer ssa, ssb
! help variables
integer nhelp1, nhelp2, nhelp4
integer nhelp41, nhelp42, nhelp21, nhelp22
integer ntest1, ntest2
integer sa1, sa2, sa3, sa4, sb1, sb2, sa34, sa134
integer ia, ib, iy, ix
integer possct

! structure A(pq,rs)*B(rs)=YC(pq)

!1.0 prepare mapdc,mapic

if ((mapda(0,6) == 1) .or. (mapda(0,6) == 4)) then
  ntest1 = 1
else
  ntest1 = 0
end if

call cct3_grc0(2,ntest1,mapda(0,1),mapda(0,2),0,0,mmul(ssa,ssb),possc0,possct,mapdc,mapic)

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
  if ((ntest2 == 1) .and. (sb1 < sb2)) then
    ! Meggie out
    goto 100
  end if

  do sa1=1,nsym
    sa134 = mmul(sa1,sa34)

    sa2 = mmul(ssa,sa134)
    if ((ntest1 == 1) .and. (sa1 < sa2)) then
      ! Meggie out
      goto 50
    end if

    !1.3 def mvec,mapdc and mapdi

    ia = mapia(sa1,sa2,sa3)
    ib = mapib(sb1,1,1)
    iy = mapic(sa1,1,1)

    ! yes/no
    if ((mapda(ia,2) > 0) .and. (mapdb(ib,2) > 0)) then
      nhelp1 = 1
    else
      goto 50
    end if

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

50  continue
  end do
100 continue
end do
ix = ix-1

return

end subroutine cct3_grc42y
