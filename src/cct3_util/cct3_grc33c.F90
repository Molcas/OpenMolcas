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

subroutine cct3_grc33C(mapda,mapdb,mapdc,mapia,mapib,mapic,mvec,ssa,ssb,pbar,posc0,ix)

use CCT3_global, only: dimm, mmul, nsym
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: mapda(0:512,6), mapdb(0:512,6), mapdc(0:512,6), mapia(8,8,8), mapib(8,8,8), mapic(8,8,8), mvec(4096,7), ssa, &
                     ssb, pbar, posc0, ix
integer(kind=iwp) :: ia, ib, ic, nhelp1, nhelp2, nhelp21, nhelp22, nhelp3, nhelp31, nhelp32, nhelp4, nhelp41, nhelp42, nsyma2, &
                     ntest1, ntest2, posct, sa1, sa12, sa2, sa3, sb1, sb12, sb2, sb3

!1*

if (pbar == 1) then

  ! structure A(p,qr)*B(qr,s)=C(p,s)

  !1.0 prepare mapdc,mapic

  call cct3_grc0(2,0,mapda(0,1),mapdb(0,3),0,0,mmul(ssa,ssb),posc0,posct,mapdc,mapic)

  !1.1 define limitations - p,q>r must be tested - ntest1

  if (mapda(0,6) == 2) then
    ntest1 = 1
  else
    ntest1 = 0
  end if

  !1.2 def symm states and test the limitations

  ix = 1
  do sa1=1,nsym

    do sa2=1,nsym
      sa12 = mmul(sa1,sa2)
      sb1 = sa2

      sa3 = mmul(ssa,sa12)
      ! Meggie out
      if ((ntest1 == 1) .and. (sa2 < sa3)) cycle
      sb2 = sa3
      sb12 = mmul(sb1,sb2)

      sb3 = mmul(ssb,sb12)

      !1.3 def mvec,mapdc and mapdi

      ia = mapia(sa1,sa2,1)
      ib = mapib(sb1,sb2,1)
      ic = mapic(sa1,sb3,1)

      ! yes/no
      if ((mapda(ia,2) > 0) .and. (mapdb(ib,2) > 0)) then
        nhelp1 = 1
      else
        cycle
      end if

      ! rowA
      nhelp2 = dimm(mapda(0,1),sa1)

      ! colB
      nhelp3 = dimm(mapdb(0,3),sb3)

      ! sum
      nhelp41 = dimm(mapda(0,2),sa2)
      nhelp42 = dimm(mapda(0,3),sa3)
      if ((ntest1 == 1) .and. (sa2 == sa3)) then
        nhelp4 = nhelp41*(nhelp41-1)/2
      else
        nhelp4 = nhelp41*nhelp42
      end if

      mvec(ix,1) = nhelp1
      mvec(ix,2) = mapda(ia,1)
      mvec(ix,3) = mapdb(ib,1)
      mvec(ix,4) = mapdc(ic,1)
      mvec(ix,5) = nhelp2
      mvec(ix,6) = nhelp4
      mvec(ix,7) = nhelp3

      ix = ix+1

    end do
  end do

else if (pbar == 2) then

  ! structure A(pq,r)*B(r,st)=C(pq,st)

  !2.1 define limitations - A p>q,r must be tested - ntest1
  !                       - B p,q>r must be tested - ntest2

  if (mapda(0,6) == 1) then
    ntest1 = 1
  else
    ntest1 = 0
  end if

  if (mapdb(0,6) == 2) then
    ntest2 = 1
  else
    ntest2 = 0
  end if

  !2.0 prepare mapdc,mapic

  if ((ntest1 == 1) .and. (ntest2 == 1)) then
    nhelp1 = 4
  else if (ntest1 == 1) then
    nhelp1 = 1
  else if (ntest2 == 1) then
    nhelp1 = 3
  else
    nhelp1 = 0
  end if

  call cct3_grc0(4,nhelp1,mapda(0,1),mapda(0,2),mapdb(0,2),mapdb(0,3),mmul(ssa,ssb),posc0,posct,mapdc,mapic)

  !2.2 def symm states and test the limitations

  ix = 1
  do sa1=1,nsym
    if (ntest1 == 1) then
      nsyma2 = sa1
    else
      nsyma2 = nsym
    end if

    do sa2=1,nsyma2
      sa12 = mmul(sa1,sa2)

      sa3 = mmul(ssa,sa12)
      sb1 = sa3

      do sb2=1,nsym
        sb12 = mmul(sb1,sb2)

        sb3 = mmul(ssb,sb12)
        ! Meggie out
        if ((ntest2 == 1) .and. (sb2 < sb3)) cycle

        !2.3 def mvec,mapdc and mapdi

        ia = mapia(sa1,sa2,sa3)
        ib = mapib(sb1,sb2,sb3)
        ic = mapic(sa1,sa2,sb2)

        ! yes/no
        if ((mapda(ia,2) > 0) .and. (mapdb(ib,2) > 0)) then
          nhelp1 = 1
        else
          cycle
        end if

        ! rowA
        nhelp21 = dimm(mapda(0,1),sa1)
        nhelp22 = dimm(mapda(0,2),sa2)
        if ((ntest1 == 1) .and. (sa1 == sa2)) then
          nhelp2 = nhelp21*(nhelp21-1)/2
        else
          nhelp2 = nhelp21*nhelp22
        end if

        ! colB
        nhelp31 = dimm(mapdb(0,2),sb2)
        nhelp32 = dimm(mapdb(0,3),sb3)
        if ((ntest2 == 1) .and. (sb2 == sb3)) then
          nhelp3 = nhelp31*(nhelp31-1)/2
        else
          nhelp3 = nhelp31*nhelp32
        end if

        ! sum
        nhelp4 = dimm(mapda(0,3),sa3)

        mvec(ix,1) = nhelp1
        mvec(ix,2) = mapda(ia,1)
        mvec(ix,3) = mapdb(ib,1)
        mvec(ix,4) = mapdc(ic,1)
        mvec(ix,5) = nhelp2
        mvec(ix,6) = nhelp4
        mvec(ix,7) = nhelp3

        ix = ix+1

      end do
    end do
  end do

end if
ix = ix-1

return

end subroutine cct3_grc33C
